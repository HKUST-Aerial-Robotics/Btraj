/*! \class GMM
    \brief Templated class which computes Group Marching Method (GMM).
    
    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to be an FMCell or something inherited from it.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
        S. Kim, An O(N) Level Set Method for Eikonal Equations, SIAM J. Sci. Comput., 22(6), 2178–2193. 2006.
        <a href="http://epubs.siam.org/doi/abs/10.1137/S1064827500367130">[PDF]</a>

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

#ifndef GMM_HPP_
#define GMM_HPP_

#include <fast_methods/fm/eikonalsolver.hpp>

template < class grid_t > class GMM : public EikonalSolver <grid_t> {

    public:
        GMM(double dt = -1) : EikonalSolver<grid_t>("GMM"), deltau_(dt) {}
        
        GMM(const char * name, double dt = -1) : EikonalSolver<grid_t>(name), deltau_(dt) {}

        virtual ~GMM() { clear(); }

        void setup
        () {
            EikonalSolver<grid_t>::setup();
            if (deltau_ < 0)
                // Original paper dt:
                //deltau_ = grid_->getLeafSize()/(grid_->getMaxSpeed()*sqrt(grid_t::getNDims()));

                // FIM paper dt:
                deltau_ = 1/grid_->getMaxSpeed();
            }

        /** \brief Actual method that implements GMM. */
        virtual void computeInternal
        () {
            if (!setup_)
                setup();

            unsigned int n_neighs;
            unsigned int j = 0;
            bool stopWavePropagation = false;

            // Algorithm initialization
            tm_= std::numeric_limits<double>::infinity();
            for (unsigned int &i: init_points_) { // For each initial point
                grid_->getCell(i).setArrivalTime(0);
                grid_->getCell(i).setState(FMState::FROZEN);
                n_neighs = grid_->getNeighbors(i, neighbors_);
                for (unsigned int s = 0; s < n_neighs; ++s){  // For each neighbor
                    j = neighbors_[s];
                    if ((grid_->getCell(j).getState() == FMState::FROZEN) || grid_->getCell(j).isOccupied())
                        continue;
                    else {
                        double new_arrival_time = solveEikonal(j);
                        if (new_arrival_time < tm_){
                            tm_ = new_arrival_time;
                        }
                        grid_->getCell(j).setArrivalTime(new_arrival_time);
                        grid_->getCell(j).setState(FMState::NARROW);
                        gamma_.push_back(j);
                    } // neighbors open.
                } // For each neighbor.
            } // For each initial point.

            // Main loop
            while(!stopWavePropagation && !gamma_.empty()) {

                tm_ += deltau_;

                std::list<unsigned int>::reverse_iterator k = gamma_.rbegin();
                std::list<unsigned int>::iterator i = k.base();//iterator points to the next element the reverse_iterator is currently pointing to
                i--;
                k = gamma_.rend();
                std::list<unsigned int>::iterator q = k.base();
                q--;//the end of a reverse list is the first element of that list
                //This is needed because some functions, like std::list::erase, do not work with reverse iterators

                // First pass
                for( ; i!=q; --i) {//for each gamma in the reverse order
                    if( grid_->getCell(*i).getArrivalTime() <= tm_) {
                        n_neighs = grid_->getNeighbors(*i, neighbors_);
                        for (unsigned int s = 0; s < n_neighs; ++s){  // For each neighbor of gamma
                            j = neighbors_[s];
                            if ( (grid_->getCell(j).getState() == FMState::FROZEN) || grid_->getCell(j).isOccupied() || (grid_->getCell(j).getVelocity() == 0)) // If Frozen,obstacle or velocity = 0
                                continue;
                            else {
                                double new_arrival_time = solveEikonal(j);
                                if (new_arrival_time < grid_->getCell(j).getArrivalTime()) // Updating narrow band if necessary.
                                    grid_->getCell(j).setArrivalTime(new_arrival_time);
                            }
                        }//for each neighbor of gamma
                    }
                }//for each gamma in the reverse order

                // Second pass
                const size_t narrow_size = gamma_.size();
                i = gamma_.begin();
                for(size_t z = 0; z < narrow_size; ++z) {//for each gamma in the forward order
                    if( grid_->getCell(*i).getArrivalTime()<= tm_) {
                        n_neighs = grid_->getNeighbors(*i, neighbors_);
                        for (unsigned int s = 0; s < n_neighs; ++s) {// for each neighbor of gamma
                            j = neighbors_[s];
                            if ((grid_->getCell(j).getState() == FMState::FROZEN) || grid_->getCell(j).isOccupied() || (grid_->getCell(j).getVelocity() == 0)) // If Frozen,obstacle or velocity = 0
                                continue;
                            else {
                                double new_arrival_time = solveEikonal(j);
                                if (new_arrival_time < grid_->getCell(j).getArrivalTime()) {
                                        grid_->getCell(j).setArrivalTime(new_arrival_time);
                                }
                                if (grid_->getCell(j).getState() == FMState::OPEN){
                                    gamma_.push_back(j);
                                    grid_->getCell(j).setState(FMState::NARROW);
                                }
                            }
                        }//for each neighbor of gamma
                    grid_->getCell(*i).setState(FMState::FROZEN);
                    if (*i == goal_idx_)
                        stopWavePropagation = true;
                    i = gamma_.erase(i);
                    }
                    else
                        ++i;
                }//for each gamma in the forward order
            }//while gamma is not zero
        }//compute

        virtual void clear
        () {
            gamma_.clear();
        }

        virtual void reset
        () {
            EikonalSolver<grid_t>::reset();
            gamma_.clear();
            tm_ = 0;
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
        /** \brief Global bound that determines the group of cells of gamma that will be updated in each step. */
        double                  tm_;

        /** \brief For each updating step, tm_ is increased by this value. */
        double                  deltau_;
        
        /** \brief List wich stores the narrow band of each iteration. */
        std::list<unsigned int> gamma_;
};

#endif /* GMM_H_*/
