/*! \class LSM
    \brief Implements Lock Sweeping Method.

    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to use an FMCell or derived.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
        S. Bak, J. McLaughlin, D. Renzi, Some Improvements for the Fast Sweeping Method,
        SIAM J. Sci. Comput., 32(5), 2853–2874. 2010.
        <a href="http://epubs.siam.org/doi/abs/10.1137/090749645">[More Info]</a>

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

#ifndef LSM_HPP_
#define LSM_HPP_

#include <fast_methods/fm/fsm.hpp>
#include <fast_methods/ndgridmap/fmcell.h>
#include <fast_methods/utils/utils.h>

/// \todo implement a more robust goal point stopping criterion.
template < class grid_t > class LSM : public FSM<grid_t> {

    public:
        LSM(unsigned maxSweeps = std::numeric_limits<unsigned>::max()) : FSM<grid_t>("LSM", maxSweeps) {}

        LSM(const char * name, unsigned maxSweeps = std::numeric_limits<unsigned>::max()) : FSM<grid_t>(name, maxSweeps) {}

        /** \brief Actual method that implements LSM. */
        virtual void computeInternal
        () {
            if (!setup_)
                setup();

            // FMState::FROZEN - locked and FMState::OPEN - unlocked.
            // The time this takes is negligible and if done in setup or
            // setEnvironment it can affect other planners run in the same
            // grid.
            for(size_t i = 0; i < grid_->size(); ++i)
                grid_->getCell(i).setState(FMState::FROZEN);

            // Initialization
            for (unsigned int i: init_points_) {
                grid_->getCell(i).setArrivalTime(0);
                unsigned int n_neighs = grid_->getNeighbors(i, neighbors_);
                for (unsigned int j = 0; j < n_neighs; ++j)
                    grid_->getCell(neighbors_[j]).setState(FMState::OPEN);
            }

            // Getting dimsizes and filling the other dimensions.
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
            FSM<grid_t>::reset();
            for(size_t i = 0; i < grid_->size(); ++i)
                grid_->getCell(i).setState(FMState::FROZEN);
        }

        virtual void printRunInfo
        () const {
            console::info("Lock Sweeping Method");
            std::cout << '\t' << name_ << '\n'
                      << '\t' << "Maximum sweeps: " << maxSweeps_ << '\n'
                      << '\t' << "Sweeps performed: " << sweeps_ << '\n'
                      << '\t' << "Elapsed time: " << time_ << " ms\n";
        }

    protected:
        /** \brief Actually executes one solving iteration of the LSM. */
        virtual void solveForIdx
        (unsigned idx) {
            if (grid_->getCell(idx).getState() == FMState::OPEN) {
                const double prevTime = grid_->getCell(idx).getArrivalTime();
                const double newTime = solveEikonal(idx);

                // Update time if better and unlock neighbors with higher time.
                if(utils::isTimeBetterThan(newTime, prevTime)) {
                    grid_->getCell(idx).setArrivalTime(newTime);
                    keepSweeping_ = true;
                    unsigned int n_neighs = grid_->getNeighbors(idx, neighbors_);
                    for (unsigned int i = 0; i < n_neighs; ++i)
                        if (utils::isTimeBetterThan(newTime, grid_->getCell(neighbors_[i]).getArrivalTime()))
                            grid_->getCell(neighbors_[i]).setState(FMState::OPEN);
                }
                // EXPERIMENTAL - Value not updated, it has converged
                else if(!isnan(newTime) && !isinf(newTime) && (idx == goal_idx_))
                    stopPropagation_ = true;

                grid_->getCell(idx).setState(FMState::FROZEN);
            }
        }

        // Inherited members from FSM.
        using FSM<grid_t>::grid_;
        using FSM<grid_t>::init_points_;
        using FSM<grid_t>::goal_idx_;
        using FSM<grid_t>::setup_;
        using FSM<grid_t>::setup;
        using FSM<grid_t>::name_;
        using FSM<grid_t>::time_;
        using FSM<grid_t>::recursiveIteration;
        using FSM<grid_t>::solveEikonal;
        using FSM<grid_t>::setSweep;
        using FSM<grid_t>::sweeps_;
        using FSM<grid_t>::maxSweeps_;
        using FSM<grid_t>::keepSweeping_;
        using FSM<grid_t>::stopPropagation_;
        using FSM<grid_t>::incs_;
        using FSM<grid_t>::inits_;
        using FSM<grid_t>::ends_;
        using FSM<grid_t>::d_;

        /** \brief Auxiliar array which stores the neighbor of each iteration of the computeFM() function. */
        std::array <unsigned int, 2*grid_t::getNDims()> neighbors_;
};

#endif /* LSM_HPP_*/
