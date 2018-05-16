/*! \class FIM
    \brief Implements Fast Iterative Method.
    
    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to use an FMCell or derived

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
        W. Jeong and R. Whitaker, A Fast Iterative Method for Eiknal Equations, SIAM J. Sci. Comput., 30(5), 2512–2534. 2008.
        <a href="http://epubs.siam.org/doi/abs/10.1137/060670298">[PDF]</a>

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

#ifndef FIM_HPP_
#define FIM_HPP_

#include <list>

#include <fast_methods/fm/eikonalsolver.hpp>
#include <fast_methods/utils/utils.h>

template < class grid_t > class FIM : public EikonalSolver<grid_t> {

    public:
        FIM(double error = 0) : EikonalSolver<grid_t>("FIM"), E_(error) {}
        FIM(const char * name, double error = 0) : EikonalSolver<grid_t>(name), E_(error) {}

        virtual ~FIM() { clear(); }

        /** \brief Actual method that implements FIM. */
        virtual void computeInternal
        () {
            if (!setup_)
                setup();

            double q =-1;
            double p =-1;
            unsigned int n_neighs = 0;
            unsigned int x_nb = 0;
            bool stopWavePropagation = 0;

            // Algorithm initialization.
            for (const unsigned int& i: init_points_) {
                grid_->getCell(i).setArrivalTime(0);
                grid_->getCell(i).setState(FMState::FROZEN);

                n_neighs = grid_->getNeighbors(i, neighbors_);
                for (unsigned int s = 0; s < n_neighs; ++s) {// For each neighbor
                    x_nb = neighbors_[s];
                    if ( (grid_->getCell(x_nb).getState() == FMState::OPEN) && !grid_->getCell(x_nb).isOccupied()) {
                        active_list_.push_back(x_nb);
                        grid_->getCell(x_nb).setState(FMState::NARROW);
                    }
                }
            }

            // Main loop.
            while(!stopWavePropagation && !active_list_.empty()) {
                for (std::list<unsigned int>::iterator x = active_list_.begin(); x!=active_list_.end(); ++x) { // for each cell of active_list
                    p = grid_->getCell(*x).getArrivalTime();
                    q = solveEikonal(*x);
                    grid_->getCell(*x).setArrivalTime(q);
                    if (fabs(p - q) <= E_) { // if the cell has converged
                        n_neighs = grid_->getNeighbors(*x, neighbors_);
                        for (unsigned int s = 0; s < n_neighs; ++s){ // For each neighbor of converged cells of active_list
                            x_nb = neighbors_[s];
                            if (grid_->getCell(x_nb).getState() != FMState::NARROW && !grid_->getCell(x_nb).isOccupied()) {
                                p = grid_->getCell(x_nb).getArrivalTime();
                                q = solveEikonal(x_nb);
                                if (utils::isTimeBetterThan(q, p)) {
                                    grid_->getCell(x_nb).setArrivalTime(q);
                                    active_list_.insert(x, x_nb);
                                    grid_->getCell(x_nb).setState(FMState::NARROW);
                                    }
                            }
                        }// For each neighbor of converged cells of active_list
                    if (*x == goal_idx_)
                        stopWavePropagation = true;
                    grid_->getCell(*x).setState(FMState::FROZEN);
                    x = active_list_.erase(x);
                    --x;
                    }// if the cell has converged
                }// for each cell of active_list
            }//while active_list is not empty
        }

        virtual void clear
        () {
            active_list_.clear();
        }

        virtual void reset
        () {
            EikonalSolver<grid_t>::reset();
            active_list_.clear();
        }

    protected:
        using EikonalSolver<grid_t>::grid_;
        using EikonalSolver<grid_t>::solveEikonal;
        using EikonalSolver<grid_t>::init_points_;
        using EikonalSolver<grid_t>::goal_idx_;
        using EikonalSolver<grid_t>::setup;
        using EikonalSolver<grid_t>::setup_;
        using EikonalSolver<grid_t>::neighbors_;

    private:
        /** \brief List wich stores the narrow band of each iteration. */
        std::list<unsigned int> active_list_;
        
        /** \brief Error threshold value that reveals if a cell has converged. */
        double E_;
};

#endif /* FIM_HPP_*/
