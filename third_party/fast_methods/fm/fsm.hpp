/*! \class FSM
    \brief Implements Fast Sweeping Method.

    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to use an FMCell or derived.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
        H. Zhao, A fast sweeping method for Eikonal equations, Math. Comp. 74 (2005), 603-627.
        <a href="http://www.ams.org/journals/mcom/2005-74-250/S0025-5718-04-01678-3/S0025-5718-04-01678-3.pdf">[PDF]</a>

    NOTE: The sweeping directions are inverted with respect to the paper to make implementation easier. And sweeping
    is implemented recursively (undetermined number of nested for loops) to achieve n-dimensional behaviour.

    Copyright (C) 2015 Javier V. Gomez
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

#ifndef FSM_HPP_
#define FSM_HPP_

#include <algorithm>

#include <fast_methods/fm/eikonalsolver.hpp>
#include <fast_methods/utils/utils.h>


/// \todo implement a more robust goal point stopping criterion.
template < class grid_t > class FSM : public EikonalSolver<grid_t> {

    public:
        FSM(unsigned maxSweeps = std::numeric_limits<unsigned>::max()) : EikonalSolver<grid_t>("FSM"),
            sweeps_(0),
            maxSweeps_(maxSweeps) {}

        FSM(const char * name, unsigned maxSweeps = std::numeric_limits<unsigned>::max()) : EikonalSolver<grid_t>(name),
            sweeps_(0),
            maxSweeps_(maxSweeps) {}

        /** \brief Sets and cleans the grid in which operations will be performed.
             Since a maximum number of dimensions is assumed, fills the rest with size 1. */
        virtual void setEnvironment
        (grid_t * g) {
            EikonalSolver<grid_t>::setEnvironment(g);
            // Filling the size of the dimensions...
            std::array<unsigned, grid_t::getNDims()> dimsize = g->getDimSizes();
            size_t ncells = 1;
            for (size_t i = 0; i < grid_t::getNDims(); ++i) {
                dimsize_[i] = dimsize[i];
                ncells *= dimsize[i];
                d_[i] = ncells;
            }

            Tvalues_.reserve(grid_t::getNDims());
        }

        /** \brief Executes EikonalSolver setup and other checks. */
        virtual void setup
        () {
            EikonalSolver<grid_t>::setup();
            initializeSweepArrays();
            if (int(goal_idx_) != -1)
                console::warning("Setting a goal point in FSM (and LSM) is experimental. It may lead to wrong results.");
        }

        /** \brief Actual method that implements FSM. */
        virtual void computeInternal
        () {
            if (!setup_)
                setup();

            // Initialization
            for (unsigned int i: init_points_) // For each initial point
                grid_->getCell(i).setArrivalTime(0);

            keepSweeping_ = true;
            stopPropagation_ = false;

            while (keepSweeping_ && !stopPropagation_ && sweeps_ < maxSweeps_) {
                keepSweeping_ = false;
                setSweep();
                ++sweeps_;
                recursiveIteration(grid_t::getNDims()-1);
            }
        }

        virtual void reset
        () {
            EikonalSolver<grid_t>::reset();
            sweeps_ = 0;
            initializeSweepArrays();
        }

        virtual void printRunInfo
        () const {
            console::info("Fast Sweeping Method");
            std::cout << '\t' << name_ << '\n'
                      << '\t' << "Maximum sweeps: " << maxSweeps_ << '\n'
                      << '\t' << "Sweeps performed: " << sweeps_ << '\n'
                      << '\t' << "Elapsed time: " << time_ << " ms\n";
        }

    protected:
        /** \brief Equivalent to nesting as many for loops as dimensions. For every most inner
         * loop iteration, solveForIdx() is called for the corresponding idx. */
        void recursiveIteration
        (size_t depth, int it = 0) {
            if (depth > 0) {
                for(int i = inits_[depth]; i != ends_[depth]; i += incs_[depth])
                    recursiveIteration(depth-1, it + i*d_[depth-1]);
            }
            else {
                for(int i = inits_[0]; i != ends_[0]; i += incs_[0])
                    if (!grid_->getCell(it+i).isOccupied())
                        solveForIdx(it+i);
            }
        }

        /** \brief Actually executes one solving iteration of the FSM. */
        virtual void solveForIdx
        (unsigned idx) {
            const double prevTime = grid_->getCell(idx).getArrivalTime();
            const double newTime = solveEikonal(idx);
            if(utils::isTimeBetterThan(newTime, prevTime)) {
                grid_->getCell(idx).setArrivalTime(newTime);
                keepSweeping_ = true;
            }
            // EXPERIMENTAL - Value not updated, it has converged
            else if(!isnan(newTime) && !isinf(newTime) && (idx == goal_idx_))
                stopPropagation_ = true;
        }

        /** \brief Set the sweep variables: initial and final indices for iterations,
             and the increment of each iteration in every dimension.

             Generates a periodical pattern for incs_ (example for 3D):
            [-1-1-1, 1-1-1, -11-1, 11-1, -1-11, 1-11, -111, 111]

             Stablishes inits_ and ends_ accordingly. */
        virtual void setSweep
        () {
            // Inspired in http://stackoverflow.com/a/17758788/2283531
            for (size_t i = 0; i < grid_t::getNDims(); ++i)
            {
                if((incs_[i] += 2) <=1)
                    break;
                else
                    incs_[i] = -1;
            }

            // Setting inits and ends.
            for (size_t i = 0; i < grid_t::getNDims(); ++i)
            {
                if (incs_[i] == 1)
                {
                    inits_[i] = 0;
                    ends_[i] = dimsize_[i];
                }
                else
                {
                    inits_[i] = dimsize_[i]-1;
                    ends_[i] = -1;
                }
            }
        }

        /** \brief Initializes the internal arrays employed. */
        virtual void initializeSweepArrays
        () {
            for (size_t i = 0; i < grid_t::getNDims(); ++i) {
                incs_[i] = 1;
                inits_[i] = 0;
                ends_[i] = 1;
            }
        }

        using EikonalSolver<grid_t>::grid_;
        using EikonalSolver<grid_t>::init_points_;
        using EikonalSolver<grid_t>::goal_idx_;
        using EikonalSolver<grid_t>::setup_;
        using EikonalSolver<grid_t>::name_;
        using EikonalSolver<grid_t>::time_;
        using EikonalSolver<grid_t>::Tvalues_;
        using EikonalSolver<grid_t>::solveEikonal;

        /** \brief Number of sweeps performed. */
        unsigned int sweeps_;

        /** \brief Number of maximum sweeps to perform. */
        unsigned maxSweeps_;

        /** \brief Flag to indicate that at least one more sweep is required. */
        bool keepSweeping_;

        /** \brief Flag to stop sweeping (used when goal point has converged). */
        bool stopPropagation_;

        /** \brief Sweep directions {-1,1} for each dimension. Extended dimensions always 1. */
        std::array<int, grid_t::getNDims()> incs_;

        /** \brief Initial indices for each dimension. Extended dimensions always 0. */
        std::array<int, grid_t::getNDims()> inits_;

        /** \brief Final indices for each dimension. Extended dimensions always 1. */
        std::array<int, grid_t::getNDims()> ends_;

        /** \brief Size of each dimension, extended to the maximum size. Extended dimensions always 1. */
        std::array<int, grid_t::getNDims()> dimsize_;

        /** \brief Auxiliar array to speed up indexing generalization: stores parcial multiplications of dimensions sizes. d_[0] = dimsize_[0];
            d_[1] = dimsize_[0]*dimsize_[1]; etc. */
        std::array<int, grid_t::getNDims()> d_;
};

#endif /* FSM_HPP_*/
