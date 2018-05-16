/*! \class SFMMStarStar
    \brief Implements the Simplified Fast Marching Method, encapsulating FMM with a priority queue.

    It uses as a main container the nDGridMap class. The nDGridMap type T
    has to be an FMCell or something inherited from it.

    The grid is assumed to be squared, that is Delta(x) = Delta(y) = leafsize_

    @par External documentation:
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

#ifndef SFMMSTAR_HPP_
#define SFMMSTAR_HPP_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <array>

#include <fast_methods/fm/sfmm.hpp>

template < class grid_t, class cell_t = FMCell>  class SFMMStar : public SFMM<grid_t, cell_t> {
    public:
        SFMMStar(HeurStrategy h = TIME) : SFMM<grid_t, cell_t>("SFMM*", h) {}
        SFMMStar(const char * name, HeurStrategy h = TIME) : SFMM<grid_t, cell_t>(name, h){}
};

#endif /* SFMMSTAR_HPP_*/
