/*! \class FM2
    \brief Implements the Fast Marching Square (FM2) planning algorithm.

    It uses as a main container the nDGridMap class. The nDGridMap template parameter
    has to be an FMCell or something inherited from it. It also uses a heap type in order
    to specify the underlying FMM.

    IMPORTANT NOTE: When running FM2 many times on the same grid it is recommended
    to completely restart the grid (erase and create or resize). See test_fm2.cpp.

    @par External documentation:
        FM2:
          A. Valero, J.V. Gómez, S. Garrido and L. Moreno,
          The Path to Efficiency: Fast Marching Method for Safer, More Efficient Mobile Robot Trajectories,
          IEEE Robotics and Automation Magazine, Vol. 20, No. 4, 2013.

    Copyright (C) 2014 Javier V. Gomez and Jose Pardeiro
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

#ifndef FM2_HPP_
#define FM2_HPP_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>
#include <limits>

#include <fast_methods/fm/fmm.hpp>
#include <fast_methods/gradientdescent/gradientdescent.hpp>

/// \todo Include support to other solvers (GMM, FIM, UFMM). It requires a better way of setting parameters.
//template < class grid_t, class solver_t = FMM<grid_t> > class FM2 : public Solver<grid_t> {

template < class grid_t, class heap_t = FMDaryHeap<FMCell> > class FM2 : public Solver<grid_t> {
    public:
    
        /** \brief Path type encapsulation. */
        typedef std::vector< std::array<double, grid_t::getNDims()> > path_t;

        /** \brief maxDistance sets the velocities map saturation distance in real units (before normalization). */
        FM2
        (double maxDistance = -1) : Solver<grid_t>("FM2"), maxDistance_(maxDistance) {
            solver_ = new FMM<grid_t, heap_t> ();
        }

        /** \brief maxDistance sets the velocities map saturation distance in real units (before normalization). */
        FM2
        (const char * name, double maxDistance = -1) : Solver<grid_t>(name), maxDistance_(maxDistance) {
            solver_ = new FMM<grid_t, heap_t> ();
        }

        virtual ~FM2 () { clear(); }

        /** \brief Sets the environment to run the solver and sets the sources for the velocities map computation. */
        virtual void setEnvironment
        (grid_t * g) {
            Solver<grid_t>::setEnvironment(g);
            grid_->getOccupiedCells(fm2_sources_);
            solver_->setEnvironment(grid_);
        }

        /** \brief Sets up the solver to check whether is ready to run. */
        virtual void setup
        () {
            Solver<grid_t>::setup();

            if(init_points_.size() > 1) {
                console::error("FM2-based solvers currently allow only 1 initial point.");
                exit(1);
            }

            if (fm2_sources_.empty()) {
                console::error("Map has no obstacles. FM2-based solver is not running.");
                exit(1);
            }
        }

        /** \brief Implements the actual FM2 method. */
        virtual void computeInternal
        () {
            if (!setup_)
                 setup();

            computeVelocitiesMap();

            // Reset time counter so that time_ returns the time of the second wave.
            start_ = std::chrono::steady_clock::now();

            // According to the theoretical basis the wave is expanded from the goal point to the initial point.
            std::vector <unsigned int> wave_init;
            wave_init.push_back(goal_idx_);
            unsigned int wave_goal = init_points_[0];

            solver_->setInitialAndGoalPoints(wave_init, wave_goal);
            solver_->compute();
            // Restore the actual grid status.
            grid_->setClean(false);
        }

        /** \brief Computes the velocities map of the FM2 algorithm. If  maxDistance_ != -1 then the map is saturated
            to the set value. It is then normalized: velocities in [0,1]. */
        void computeVelocitiesMap
        () {
            // Forces not to clean the grid.
            grid_->setClean(true);
            solver_->setInitialPoints(fm2_sources_);
            solver_->compute();
            time_vels_ = solver_->getTime();
            start_ = std::chrono::steady_clock::now();
            // Rescaling and saturating to relative velocities: [0,1]
            double maxValue = grid_->getMaxValue();
            double maxVelocity = 0;

            if (maxDistance_ != -1)
                maxVelocity = maxDistance_ / grid_->getLeafSize();

            for (unsigned int i = 0; i < grid_->size(); ++i) {
                double vel = grid_->getCell(i).getValue() / maxValue;

                if (maxDistance_ != -1)
                    if (vel < maxVelocity)
                        grid_->getCell(i).setVelocity(vel / maxVelocity);
                    else
                        grid_->getCell(i).setVelocity(1);
                else
                    grid_->getCell(i).setVelocity(vel);

                // Restarting grid values for second wave expasion.
                grid_->getCell(i).setValue(std::numeric_limits<double>::infinity());
                grid_->getCell(i).setState(FMState::OPEN);
                grid_->setClean(true);
            }
            end_ = std::chrono::steady_clock::now();
            time_vels_ += std::chrono::duration_cast<std::chrono::milliseconds>(end_-start_).count();
        }

        /** \brief Encapsulates the path extraction.

            Computes the path from the given goal index to the minimum
            of the times of arrival map. No checks are done (points in the borders, points in obstacles...).
            @param p path the resulting path (output).
            @param path_velocity velocity the resulting path (output).
            @param goal_idx index of the goal point, where gradient descent will start. If
                   no specified, the previously set goal point is used. */
        virtual void computePath
        (path_t * p, std::vector <double> * path_velocity, double step = 1) {
            path_t* path_ = p;
            GradientDescent< nDGridMap<FMCell, grid_t::getNDims()> > grad;
            grad.apply(*grid_,init_points_[0],*path_, *path_velocity, step);
        }

        virtual void clear
        () {
            Solver<grid_t>::clear();
            fm2_sources_.clear();
            maxDistance_ = -1;
            delete solver_;
        }

        virtual void reset
        () {
            Solver<grid_t>::reset();
            solver_->reset();
        }

        /** \brief Returns velocities map computation time. */
        virtual double getTimeVelocities
        () const {
            return time_vels_;
        }

    protected:
        using Solver<grid_t>::grid_;
        using Solver<grid_t>::init_points_;
        using Solver<grid_t>::goal_idx_;
        using Solver<grid_t>::setup_;
        using Solver<grid_t>::time_;
        using Solver<grid_t>::start_;
        using Solver<grid_t>::end_;

        /** \brief Wave propagation sources for the Fast Marching Square velocities map computation.*/
        std::vector<unsigned int>   fm2_sources_;
        
        /** \brief Underlying FMM-based solver. */
        FMM<grid_t, heap_t> *       solver_;
        
        /** \brief Distance value to saturate the first potential. */
        double                      maxDistance_;

        /** \brief Time elapsed in the velocities map computation. */
        double                      time_vels_;
};

#endif /* FM2_H_*/
