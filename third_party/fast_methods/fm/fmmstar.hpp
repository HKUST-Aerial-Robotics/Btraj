/*! \class FMMStar
    \brief Encapsulates the calls to FMM with heuristics enabled.

    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to be an FMCell or something inherited from it.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    The type of the heap introduced is very important for the behaviour of the
    algorithm. The following heaps are provided:

    - FMDaryHeap wrap for the Boost D_ary heap (generalization of binary heaps).
    * Set by default if no other heap is specified. The arity has been set to 2
    * (binary heap) since it has been tested to be the more efficient in this algorithm.
    - FMFibHeap wrap for the Boost Fibonacci heap.
    - FMPriorityQueue wrap to the std::PriorityQueue class. This heap implies the implementation
    * of the Simplified FMMStar (SFMMStar) method, done automatically because of the FMPriorityQueue::increase implementation.
    *
    @par External documentation:
        FMMStar:
          A. Valero, J.V. Gómez, S. Garrido and L. Moreno, The Path to Efficiency: Fast Marching Method for Safer, More Efficient Mobile Robot Trajectories, IEEE Robotics and Automation Magazine, Vol. 20, No. 4, 2013. DOI: <a href="http://dx.doi.org/10.1109/MRA.2013.2248309">10.1109/MRA.2013.2248309></a><br>
           <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6582543">[PDF]</a>

        SFMMStar:
          M.W. Jones, J.A. Baerentzen, M. Sramek, 3D Distance Fields: A Survey of Techniques and Applications, IEEE Transactions on Visualization and Computer Graphics, Vol. 12, No. 4, 2006. DOI <a href=http://dx.doi.org/10.1109/TVCG.2006.56">110.1109/TVCG.2006.56</a><br>
          <a href="http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=1634323">[PDF]</a>

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

#ifndef FMMStar_HPP_
#define FMMStar_HPP_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>

#include <fast_methods/fm/fmm.hpp>

#include <fast_methods/ndgridmap/fmcell.h>
#include <fast_methods/datastructures/fmdaryheap.hpp>

#include <fast_methods/ndgridmap/ndgridmap.hpp>
#include <fast_methods/console/console.h>

template < class grid_t, class heap_t = FMDaryHeap<FMCell> >  class FMMStar : public FMM<grid_t, heap_t> {

    /** \brief Shorthand for base solver. */
    typedef FMM<grid_t, heap_t> FMMBase;

    public:
        FMMStar(HeurStrategy h = TIME) : FMMBase("FMM*", h) {
            /// \todo automate the naming depending on the heap.
            //if (static_cast<FMFibHeap>(heap_t))
             //   name_ = "FMMStarFib";
        }

        FMMStar(const char * name, HeurStrategy h = TIME) : FMMBase(name, h){}
};

#endif /* FMMSTAR_HPP_*/
