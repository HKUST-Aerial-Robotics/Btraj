/*! \file FMDirectionalcell.h
    \brief Header of the FMDirectionalCell class
    
    A stand-alone, standard C++ class which represents each one of the cells
    of a gridmap and its typical members. Inherited from Cell class, in this
    case the value_ member represents the distance value (or time of arrival).

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

#ifndef FMDIRECTIONALCELL_H_
#define FMDIRECTIONALCELL_H_

#include <iostream>
#include <string>
#include <limits>

#include "../../ndgridmap/cell.h"
#include "fmcell.h"

class FMDirectionalCell : public FMCell{
      /**
       * ostream operator << overloaded for this class.
       */
    friend std::ostream& operator << (std::ostream & os, const FMDirectionalCell & c);

    public:
        /**
         * Default constructor which performs and implicit Fast Marching-like initialization of the grid,
         */
        FMDirectionalCell() : FMCell() {}

        virtual ~FMDirectionalCell() {}
        
        // NOTE: no checks are done (out of bounds, correct states, etc) no improve efficienty.
        // TODO: overload functions to add the option of input checking.
        virtual void setDirectionalTime (const double dir_time) { dir_time_= dir_time; }

        virtual double getDirectionalTime () const                {return dir_time_;}

    protected:

        double dir_time_;
};

#endif /* FMDIRECTIONALCELL_H_*/
