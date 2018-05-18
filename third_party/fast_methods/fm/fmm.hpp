/*! \class FMM
    \brief Implements the Fast Marching Method (FMM).

    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to be an FMCell or something inherited from it.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    The type of the heap introduced is very important for the behaviour of the
    algorithm. The following heaps are provided:

    - FMDaryHeap wrap for the Boost D_ary heap (generalization of binary heaps).
    * Set by default if no other heap is specified. The arity has been set to 2
    * (binary heap) since it has been tested to be the more efficient in this algorithm.
    - FMFibHeap wrap for the Boost Fibonacci heap.
    - FMPriorityQueue wrap to the std::PriorityQueue class. This heap implies the implementation
    * of the Simplified FMM (SFMM) method, done automatically because of the FMPriorityQueue::increase implementation.

    @par External documentation:
        FMM:
          A. Valero, J.V. Gómez, S. Garrido and L. Moreno, The Path to Efficiency: Fast Marching Method for Safer, More Efficient Mobile Robot Trajectories, IEEE Robotics and Automation Magazine, Vol. 20, No. 4, 2013. DOI: <a href="http://dx.doi.org/10.1109/MRA.2013.2248309">10.1109/MRA.2013.2248309></a><br>
           <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6582543">[PDF]</a>

        SFMM:
          M.W. Jones, J.A. Baerentzen, M. Sramek, 3D Distance Fields: A Survey of Techniques and Applications, IEEE Transactions on Visualization and Computer Graphics, Vol. 12, No. 4, 2006. DOI <a href=http://dx.doi.org/10.1109/TVCG.2006.56">110.1109/TVCG.2006.56</a><br>
          <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=1634323">[PDF]</a>

    Copyright (C) 2014 Javier V. Gomez
    www.javiervgomez.com

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FMM_HPP_
#define FMM_HPP_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>

#include <fast_methods/fm/eikonalsolver.hpp>

#include <fast_methods/ndgridmap/fmcell.h>
#include <fast_methods/datastructures/fmdaryheap.hpp>

#include <fast_methods/ndgridmap/ndgridmap.hpp>
#include <fast_methods/console/console.h>

/** \brief Heuristic strategy to be used. TIME = DISTANCE/local velocity. */
enum HeurStrategy {NOHEUR = 0, TIME, DISTANCE};

template < class grid_t, class heap_t = FMDaryHeap<FMCell> >  class FMM : public EikonalSolver<grid_t> {

    public:
        FMM(HeurStrategy h = NOHEUR) : EikonalSolver<grid_t>("FMM"), heurStrategy_(h), precomputed_(false) {
            /// \todo automate the naming depending on the heap.
            //if (static_cast<FMFibHeap>(heap_t))
             //   name_ = "FMMFib";
        }

        FMM(const char * name, HeurStrategy h = NOHEUR) : EikonalSolver<grid_t>(name), heurStrategy_(h), precomputed_(false) {}

        virtual ~FMM() { clear(); }

        /** \brief Executes EikonalSolver setup and sets maximum size for the narrow band. */
        virtual int setup
        () {
            int ret = EikonalSolver<grid_t>::setup();
            narrow_band_.setMaxSize(grid_->size());
            setHeuristics(heurStrategy_); // Redundant, but safe.

            if (int(goal_idx_) == -1 && heurStrategy_ != NOHEUR) {
                console::warning("FMM/SFMM: Heuristics set with no goal point. Deactivating heuristics.");
                heurStrategy_ = NOHEUR;
            }

            return ret;
        }

        /** \brief Actual method that implements FMM. */
        virtual int computeInternal(double max_v) 
        {
            if (!setup_)
                if(setup() == -1)
                    return -1;
            
            unsigned int j = 0;
            unsigned int n_neighs = 0;
            bool stopWavePropagation = false;

            // Algorithm initialization
            for (unsigned int &i: init_points_) 
            { // For each initial point
                grid_->getCell(i).setArrivalTime(0);
                // Include heuristics if necessary.
                if (heurStrategy_ == TIME)
                {   
                    grid_->getCell(i).setHeuristicTime( getPrecomputedDistance(i) / max_v);//grid_->getCell(i).getVelocity() );
                }
                else if (heurStrategy_ == DISTANCE)
                    grid_->getCell(i).setHeuristicTime( getPrecomputedDistance(i) );
                narrow_band_.push( &(grid_->getCell(i)) );
            }

            // Main loop.
            unsigned int idxMin = 0;
            int iter = 0;
            while (!stopWavePropagation && !narrow_band_.empty()) 
            {   
                iter++;
                idxMin = narrow_band_.popMinIdx();

                n_neighs = grid_->getNeighbors(idxMin, neighbors_);
                grid_->getCell(idxMin).setState(FMState::FROZEN);
                for (unsigned int s = 0; s < n_neighs; ++s) 
                {
                    j = neighbors_[s];

                    if ((grid_->getCell(j).getState() == FMState::FROZEN) || grid_->getCell(j).isOccupied())
                        continue;
                    else 
                    {
                        double new_arrival_time = solveEikonal(j);

                        // Include heuristics if necessary.
                        if (heurStrategy_ == TIME)
                            grid_->getCell(j).setHeuristicTime( getPrecomputedDistance(j) / max_v);//grid_->getCell(j).getVelocity() );
                        else if (heurStrategy_ == DISTANCE)
                            grid_->getCell(j).setHeuristicTime( getPrecomputedDistance(j) );

                        // Updating narrow band if necessary.
                        if (grid_->getCell(j).getState() == FMState::NARROW) 
                        {
                            if (utils::isTimeBetterThan(new_arrival_time, grid_->getCell(j).getArrivalTime())) 
                            {
                                grid_->getCell(j).setArrivalTime(new_arrival_time);
                                narrow_band_.increase( &(grid_->getCell(j)) );
                            }
                        }
                        else 
                        {
                            grid_->getCell(j).setState(FMState::NARROW);
                            grid_->getCell(j).setArrivalTime(new_arrival_time);
                            narrow_band_.push( &(grid_->getCell(j)) );
                        } // neighbors_ open.
                    } // neighbors_ not frozen.
                } // For each neighbor.

                if (idxMin == goal_idx_)
                    stopWavePropagation = true;
            } // while narrow band not empty

            //std::cout<<"iteration num: "<<iter<<endl;
            return 1;
        }

        /** \brief Set heuristics flag. True is activated. It will precompute distances
            if not done already. */
        void setHeuristics(HeurStrategy h) 
        {
            if (h && int(goal_idx_)!=-1) 
            {
                heurStrategy_ = h;
                grid_->idx2coord(goal_idx_, heur_coord_);
                //grid_->idx2coord(goal_idx_, heur_coord_);
                if (!precomputed_)
                    precomputeDistances();
            }
        }

        /** \brief Returns heuristics flag. */
        HeurStrategy getHeuristics
        () const {
            return heurStrategy_;
        }

        virtual void clear
        () {
            narrow_band_.clear();
            distances_.clear();
            precomputed_ = false;
        }

        virtual void reset
        () {
            EikonalSolver<grid_t>::reset();
            narrow_band_.clear();
        }

        /** \brief Computes euclidean distance between goal and rest of cells. */
        virtual void precomputeDistances() 
        {   
            distances_.reserve(grid_->size());
            std::array <unsigned int, grid_t::getNDims()> coords;
            double dist = 0;

            for (size_t i = 0; i < grid_->size(); ++i)
            {
                dist = 0;
                grid_->idx2coord(i, coords);

                for (size_t j = 0; j < coords.size(); ++j)
                    dist += ((int)coords[j] - (int)heur_coord_[j]) * ((int)coords[j] - (int)heur_coord_[j]);

                distances_[i] = 1.00001 * std::sqrt(dist) * grid_->getLeafSize();
                
                /*for (size_t j = 0; j < coords.size(); ++j)
                    dist += fabs((int)coords[j] - (int)heur_coord_[j]);
                distances_[i] = 1.00001 * grid_->getLeafSize();*/
                
            }
            precomputed_ = true;
        }

        /** \brief Extracts the euclidean distance calculated from precomputeDistances
            function distance between two positions. */
        virtual double getPrecomputedDistance(const unsigned int idx) 
        {
            return distances_[idx];
        }

        virtual void printRunInfo
        () const {
            console::info("Fast Marching Method");
            std::cout << '\t' << name_ << '\n'
                      << '\t' << "Heuristic type: " << heurStrategy_ << '\n'
                      << '\t' << "Elapsed time: " << time_ << " ms\n";
        }


    /// \note These accessing levels may need to be modified (and other EikonalSolvers).
    protected:
        using EikonalSolver<grid_t>::grid_;
        using EikonalSolver<grid_t>::init_points_;
        using EikonalSolver<grid_t>::goal_idx_;
        using EikonalSolver<grid_t>::setup_;
        using EikonalSolver<grid_t>::name_;
        using EikonalSolver<grid_t>::time_;
        using EikonalSolver<grid_t>::solveEikonal;
        using EikonalSolver<grid_t>::neighbors_;

    private:
        /** \brief Instance of the heap used. */
        heap_t                                          narrow_band_;

        /** \brief Flag to activate heuristics and corresponding strategy. */
        HeurStrategy                                    heurStrategy_;

        /** \brief Stores the precomputed heuristic distances. */
        std::vector<double>                             distances_;
        
        /** \brief Flag to indicate if distances_ is already computed. */
        bool                                            precomputed_;

        /** \brief Goal coord, goal of the second wave propagation (actually the initial point of the path). */
        std::array <unsigned int, grid_t::getNDims()>   heur_coord_;
};

#endif /* FMM_HPP_*/
