/*! \class FMCell
    \brief A stand-alone, standard C++ class which represents each one of the cells
    of a gridmap and its typical members when used together with Fast Marching Methods.
    Inherited from Cell class, in this case the value_ member represents the distance value (or time of arrival)
    and occupancy_ represents the propagation velocity. 
   
    IMPORTANT NOTE: no checks are done in the set functions.
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

#ifndef FMCELL_H_
#define FMCELL_H_

#include <iostream>
#include <string>
#include <limits>

#include "../../ndgridmap/cell.h"

/** \brief Possible states of the FMCells*/
enum class FMState {OPEN, NARROW, FROZEN};

/// \todo Overload functions to add the option of input checking. No checks are faster.
class FMCell : public Cell{
    friend std::ostream& operator << (std::ostream & os, const FMCell & c);

    public:
        /** \brief Default constructor which performs and implicit Fast Marching-like initialization of the grid. */
        FMCell() : Cell(std::numeric_limits<double>::infinity(), 1), state_(FMState::OPEN), bucket_(0), hValue_(0) {}

        virtual ~FMCell() {}

        virtual inline void setVelocity(double v)           {occupancy_ = v;}
        virtual inline void setArrivalTime(double at)       {value_= at;}
        virtual inline void setHeuristicTime(double hv)     {hValue_ = hv; }//std::cout<<"h: "<<hv<<std::endl; }
        virtual inline void setState(FMState state)         {state_ = state;}
        virtual inline void setBucket(int b)                {bucket_ = b;}
        
        /** \brief Sets default values for the cell. Concretely, restarts value_ = Inf, state_ = OPEN and
            hValue_ = 0 but occupancy_ is not modified. */
        virtual void setDefault();

        std::string type() {return std::string("FMCell - Fast Marching cell");}

        virtual inline double getArrivalTime() const              {return value_;}
        virtual inline double getHeuristicValue() const           {return hValue_;}
        virtual inline double getTotalValue() const               
        {
            return value_ + 2.0 * hValue_; // in very large scale with high obs ddensity, h > 1 improve performance sigfcantly
        }
        virtual inline double getVelocity() const                 {return occupancy_;}
        virtual inline FMState getState() const                   {return state_;}
        virtual inline int getBucket() const                      {return bucket_;}

    protected:
        /** \brief State of the cell. */
        FMState state_;

        /** \brief Used when sorted with FMUntidyQueue. */
        int bucket_;

        /** \brief Heuristic value. */
        double hValue_;
};

#endif /* FMCELL_H_*/