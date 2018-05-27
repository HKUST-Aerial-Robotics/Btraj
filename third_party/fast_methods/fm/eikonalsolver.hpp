/*! \class EikonalSolver
    \brief Abstract class that serves as interface for the actual EikonalSolvers implemented.
    It requires (at least) the computeInternal method to be implemented,

    It uses as a main container the nDGridMap class. The nDGridMap template paramenter
    has to be an FMCell or something inherited from it.

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

#ifndef EIKONALSOLVER_H_
#define EIKONALSOLVER_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>
#include <chrono>

#include <boost/concept_check.hpp>

#include <fast_methods/fm/solver.hpp>
#include <fast_methods/console/console.h>

template <class grid_t>
class EikonalSolver : public Solver<grid_t>{

    public:
        EikonalSolver() : Solver<grid_t>("EikonalSolver") {}
        EikonalSolver(const std::string& name) : Solver<grid_t>(name) {}

        /** \brief Solves nD Eikonal equation for cell idx. If heuristics are activated, it will add
            the estimated travel time to goal with current velocity. */
        virtual double solveEikonal(const int & idx) 
        {   
            unsigned int a = grid_t::getNDims(); // a parameter of the Eikonal equation.
            Tvalues_.clear();

            for (unsigned int dim = 0; dim < grid_t::getNDims(); ++dim) {
                double minTInDim = grid_->getMinValueInDim(idx, dim);
                if (!std::isinf(minTInDim) && minTInDim < grid_->getCell(idx).getArrivalTime())
                    Tvalues_.push_back(minTInDim);
                else
                    a -=1;
            }

            if (a == 0)
                return std::numeric_limits<double>::infinity();

            // Sort the neighbor values to make easy the following code.
            /// \todo given that this sorts a small vector, a n^2 methods could be better. Test it.
            std::sort(Tvalues_.begin(), Tvalues_.end());
            double updatedT;
            for (unsigned i = 1; i <= a; ++i) {
                updatedT = solveEikonalNDims(idx, i);
                // If no more dimensions or increasing one dimension will not improve time.
                if (i == a || (updatedT - Tvalues_[i]) < utils::COMP_MARGIN)
                    break;
            }
            return updatedT;
        }

    protected:
        /** \brief Solves the Eikonal equation assuming that Tvalues_
            is sorted. */
        double solveEikonalNDims
        (unsigned int idx, unsigned int dim) {
            // Solve for 1 dimension.
            if (dim == 1)
                return Tvalues_[0] + grid_->getLeafSize() / grid_->getCell(idx).getVelocity();

            // Solve for any number > 1 of dimensions.
            double sumT = 0;
            double sumTT = 0;
            for (unsigned i = 0; i < dim; ++i) {
                sumT += Tvalues_[i];
                sumTT += Tvalues_[i]*Tvalues_[i];
            }

            // These a,b,c values are simplified since leafsize^2, which should be present in the three
            // terms but they are cancelled out when solving the quadratic function.
            double a = dim;
            double b = -2*sumT;
            double c = sumTT - grid_->getLeafSize() * grid_->getLeafSize() / (grid_->getCell(idx).getVelocity()*grid_->getCell(idx).getVelocity());
            double quad_term = b*b - 4*a*c;

            if (quad_term < 0)
                return std::numeric_limits<double>::infinity();
            else
                return (-b + sqrt(quad_term))/(2*a);
        }

        /** \brief Auxiliar vector with values T0,T1...Tn-1 variables in the Discretized Eikonal Equation. */
        std::vector<double>          Tvalues_;

        /** \brief Auxiliar array which stores the neighbor of each iteration of the computeFM() function. */
        std::array <unsigned int, 2*grid_t::getNDims()> neighbors_;

        using Solver<grid_t>::grid_;
};

#endif /* EIKONALSOLVER_H_*/
