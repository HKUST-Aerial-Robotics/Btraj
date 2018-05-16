/*! \class UFMM
    \brief Fast Marching Method using a untidy priority queue (UFMM).
    
    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to be an FMUntidyCell or something inherited from it.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
        L. Yatziv, A.Bartesaghi and G. Sapiro, O(n) implementation of the fast marching algorithm, Journal of Computational Physics. 212(2): 393-399. 2006.
        <a href="http://www.sciencedirect.com/science/article/pii/S0021999105003736">[PDF]</a>
         
        C. Rasch, T. Satzger, Remarks on the O(N) implementation of the Fast Marching Method. 2007.
        <a href="http://arxiv.org/abs/cs/0703082">[PDF]</a>
    
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

#ifndef UFMM_HPP_
#define UFMM_HPP_

#include <fast_methods/fm/eikonalsolver.hpp>
#include <fast_methods/datastructures/fmuntidyqueue.hpp>

template <class grid_t, class cell_t = FMCell> class UFMM : public EikonalSolver<grid_t> {

    public:
        UFMM
        (unsigned s = 1000, double inc = 2) : EikonalSolver<grid_t>("UFMM"), heap_s_(s), heap_inc_(inc) {
            narrow_band_ = new FMUntidyQueue<cell_t> (heap_s_, heap_inc_);
        }
        UFMM
        (const char * name, unsigned s = 1000, double inc = 2) : EikonalSolver<grid_t>(name), heap_s_(s), heap_inc_(inc) {
            narrow_band_ = new FMUntidyQueue<cell_t> (heap_s_, heap_inc_);
        }

        virtual ~UFMM() { clear(); }

        /** \brief Actual method that implements UFMM. */
        virtual void computeInternal
        () {
            if (!setup_)
                setup();

            unsigned int j= 0;
            unsigned int n_neighs = 0;
            bool stopWavePropagation = false;

            // Algorithm initialization
            for (unsigned int &i : init_points_) { // For each initial point
                grid_->getCell(i).setArrivalTime(0);
                narrow_band_->push( &(grid_->getCell(i)) );
            }

            // Main loop.
            unsigned int idxMin = 0;
            while (!stopWavePropagation && !narrow_band_->empty()) {
                idxMin = narrow_band_->topIdx(); // pop() has to be called after pushing in this case (because
                                                 // of the untidy queue implementation.
                grid_->getCell(idxMin).setState(FMState::FROZEN);
                n_neighs = grid_->getNeighbors(idxMin, neighbors_);
                for (unsigned int s = 0; s < n_neighs; ++s) { // For each neighbor.
                    j = neighbors_[s];
                    if ( (grid_->getCell(j).getState() == FMState::FROZEN) || grid_->getCell(j).isOccupied())
                        continue;
                    else {
                        double new_arrival_time = solveEikonal(j);
                        if (grid_->getCell(j).getState() == FMState::NARROW) { // Updating narrow band if necessary.
                            if (utils::isTimeBetterThan(new_arrival_time, grid_->getCell(j).getArrivalTime()) ) {
                                grid_->getCell(j).setArrivalTime(new_arrival_time);
                                narrow_band_->increase( &(grid_->getCell(j)) );
                            }
                        }
                        else {
                            grid_->getCell(j).setState(FMState::NARROW);
                            grid_->getCell(j).setArrivalTime(new_arrival_time);
                            narrow_band_->push( &(grid_->getCell(j)) );
                        } // neighbors open.
                    } // neighbors not frozen.
                } // For each neighbor.
                narrow_band_->pop();
                if (idxMin == goal_idx_)
                    stopWavePropagation = true;
            } // while narrow band is not empty
        }

        virtual void printRunInfo
        () const {
            console::info("Untidy Fast Marching Method");
            std::cout << '\t' << name_ << '\n'
                      << '\t' << "Number of buckets: " << heap_s_ << '\n'
                      << '\t' << "Maximum increment " << heap_inc_ << '\n'
                      << '\t' << "Elapsed time: " << time_ << " ms\n";
        }

        virtual void clear
        () {
            delete narrow_band_;
        }

        virtual void reset
        () {
            EikonalSolver<grid_t>::reset();
            narrow_band_->clear();
        }

    protected:
        using EikonalSolver<grid_t>::grid_;
        using EikonalSolver<grid_t>::solveEikonal;
        using EikonalSolver<grid_t>::init_points_;
        using EikonalSolver<grid_t>::goal_idx_;
        using EikonalSolver<grid_t>::setup;
        using EikonalSolver<grid_t>::setup_;
        using EikonalSolver<grid_t>::neighbors_;
        using EikonalSolver<grid_t>::name_;
        using EikonalSolver<grid_t>::time_;

    private:
        /** \brief Number of buckets in the heap. */
        unsigned                heap_s_;

        /** \brief Size (maximum increment) of the heap. */
        double                  heap_inc_;

        /** \brief Heap Instance of the priority queue used. */
        FMUntidyQueue<cell_t> * narrow_band_;
};

#endif /* UFMM_HPP_*/
