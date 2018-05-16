/*! \struct FMCompare
    \brief Compare function to be used with FM-based heaps, concretely
    FMDaryHeap, FMFibHeap and FMPriorityQueue.
    
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

#ifndef FMCOMPARE_H_
#define FMCOMPARE_H_

#include "fmcell.h"

/** \brief This struct is used a comparator for the heap. Since a minimum-heap
    is desired the operation checked is param1 > param2 as seen in this
    [Stack Overflow post](http://stackoverflow.com/a/16706002/2283531). */
template <class cell_t> struct FMCompare {
    inline bool operator()
    (const cell_t * c1 , const cell_t * c2) const {
        return c1->getTotalValue() > c2->getTotalValue();
    }
};

#endif /* FMCOMPARE_H_ */
