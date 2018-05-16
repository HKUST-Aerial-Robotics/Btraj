/*! \class GridPlotter
    \brief Auxiliar class which helps to visualise Fast Marching steps and results.
    
    It is based on the CImg library, therefore it has to be accessible.
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

#ifndef GRIDPLOTTER_H_
#define GRIDPLOTTER_H_

#include <string>
#include <array>

#include "../fmm/fmdata/fmcell.h"
#include "../ndgridmap/ndgridmap.hpp"

#include "../thirdparty/CImg.h"

using namespace cimg_library;

/// \todo include checks which ensure that the grids are adecuate for the functions used.
/// \todo generalize plotMapPaths to work with multiple (n>2) paths.
class GridPlotter {
    /** \brief Shorthand for 2D coordinates. */
    typedef typename std::array<unsigned int, 2> Coord2D;
    
    /** \brief Shorthand for 2D real points. */
    typedef typename std::array<double, 2> Point2D;
    
    /** \brief Shorthand for 2D paths of real points. */
    typedef typename std::vector <Point2D> Path2D;
    
    /** \brief Shorthand for vector of 2D paths of real points. */
    typedef typename std::vector <Path2D> Paths2D;

    public:
        /** \brief Plots the binary map of a given grid. It is based on the
            nDGridMap::getOccupancy(). This function has to be overloaded 
            in another occupancy type is being used.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
        template<class T, size_t ndims> 
        static void plotMap
        (nDGridMap<T, ndims> & grid, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            CImg<bool> img(d[0],d[1],1,1,0);
            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXY(img,x,y) {img(x,y) = !grid[img.width()*(img.height()-y-1)+x].isOccupied(); }
            // Code for not-inverted Y axis.
            //cimg_forXY(img,x,y) { if(grid[img.width()*y+x].getOccupancy() > 0.5) img(x,y) = true;
            //                      else img(x,y) = false; }
            name += " Map";
            img.display(name.c_str(), false);
        }

        /** \brief Plots the occupancy map of a given grid. It is based on the
            nDGridMap::getOccupancy(), This function has to be overloaded if 
            another occupancy type is being used.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
        template<class T, size_t ndims>
        static void plotOccupancyMap
        (nDGridMap<T, ndims> & grid, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            CImg<double> img(d[0],d[1],1,1,0);
            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXY(img,x,y) { img(x,y) = grid[img.width()*(img.height()-y-1)+x].getOccupancy()*255; }
            name += " Occupancy Map";
            img.display(name.c_str(), false);
        }

       /** \brief Plots the values map of a given grid. It is based on the
           nDGridMap::getValue(). This function has to be overloaded in 
           another value type is being used.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

           IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getValue() method. */
        template<class T, size_t ndims = 2>
        static void plotArrivalTimes
        (nDGridMap<T, ndims> & grid, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            double max_val = grid.getMaxValue();
            CImg<double> img(d[0],d[1],1,1,0);
            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXY(img,x,y) { img(x,y) = grid[img.width()*(img.height()-y-1)+x].getValue()/max_val*255; }
            img.map( CImg<double>::jet_LUT256() );
            name += " Grid values";
            img.display(name.c_str(), false);
        }

       /** \brief Plots the values map of a given grid. It is based on the
            nDGridMap::getValue(). This function has to be overloaded in 
            another value type is being used. Also plots a given path.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
        template<class T, size_t ndims = 2>
        static void plotMapPath
        (nDGridMap<T, ndims> & grid, const Path2D & path, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            CImg<double> img(d[0],d[1],1,3,0);

            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXYZC(img,x,y,z,c) { img(x,y,z,c) = (!grid[img.width()*(img.height()-y-1)+x].isOccupied())*255; }

            for (unsigned int i = 0; i< path.size(); ++i)
            {
                img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, 1) = 0;
                img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, 2) = 0;
            }

            name += " Map and Path";
            img.display(name.c_str(), false);
        }

       /** \brief Plots the values map of a given grid. It is based on the
            nDGridMap::getValue(). This function has to be overloaded in 
            another value type is being used.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
        template<class T, size_t ndims = 2>
        static void plotOccupancyPath
        (nDGridMap<T, ndims> & grid, const Path2D & path, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            CImg<double> img(d[0],d[1],1,3,0);
            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXYZC(img,x,y,z,c) { img(x,y,z,c) = grid[img.width()*(img.height()-y-1)+x].getOccupancy()*255; }

            for (unsigned int i = 0; i< path.size(); ++i)
            {
                img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, 1) = 0;
                img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, 2) = 0;
            }
            name += " Map and Path";
            img.display(name.c_str(), false);
        }

       /** \brief Plots the values map of a given grid. It is based on the
            nDGridMap::getValue(). This function has to be overloaded in 
            another value type is being used. Also plots 2 given paths.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
        template<class T, size_t ndims = 2>
        static void plotMapPaths
        (nDGridMap<T, ndims> & grid, const Paths2D & paths, std::string name = "") {
            std::array<unsigned int,2> d = grid.getDimSizes();
            CImg<double> img(d[0],d[1],1,3,0);

            // Filling the image flipping Y dim. We want now top left to be the (0,0).
            cimg_forXYZC(img,x,y,z,c) { img(x,y,z,c) = (!grid[img.width()*(img.height()-y-1)+x].isOccupied())*255; }

            // Draw the path using different colours
            for (unsigned int j = 0; j < paths.size(); ++j)
            {
                Path2D path = paths[j];
                for (unsigned int i = 0; i< path.size(); ++i)
                {
                    img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, j) = 0;
                    img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1), 0, j+1) = 0;
                }
            }
            name += " Map and Paths";
            img.display(name.c_str(), false);
        }

       /** \brief Plots the values map of a given grid. It is based on the
            nDGridMap::getValue(). This function has to be overloaded in 
            another value type is being used. It also plots a given path.

            Should be used only in 2D grids.

            The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

            IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getOccupancy() method. */
      template<class T, size_t ndims = 2>
      static void plotArrivalTimesPath
      (nDGridMap<T, ndims> & grid, const Path2D & path, std::string name = "") {
          std::array<unsigned int,2> d = grid.getDimSizes();
          double max_val = grid.getMaxValue();
          CImg<double> img(d[0],d[1],1,1,0);

          // Filling the image flipping Y dim. We want now top left to be the (0,0).
          cimg_forXY(img,x,y) { img(x,y) = grid[img.width()*(img.height()-y-1)+x].getValue()/max_val*255; }

          for (unsigned int i = 0; i< path.size(); ++i)
              img(static_cast<unsigned int>(path[i][0]), (img.height()-static_cast<unsigned int>(path[i][1])-1)) = 255;

          img.map( CImg<double>::jet_LUT256() );
          name += " Values and Path";
          img.display(name.c_str(), false);
      }

      /** \brief Plots the FMState of a cell.
           Should be used only in 2D grids.

           The Y dimension flipping is because nDGridMap works in X-Y coordinates, not in image indices as CImg.

          IMPORTANT NOTE: no type-checkings are done. T type has to be Cell or any class with bool getValue() method. */
       template<class T, size_t ndims = 2>
       static void plotFMStates
       (nDGridMap<T, ndims> & grid, std::string name = "") {
           std::array<unsigned int,2> d = grid.getDimSizes();
           //double max_val = grid.getMaxValue();
           CImg<unsigned int> img(d[0],d[1],1,1,0);
           // Filling the image flipping Y dim. We want now top left to be the (0,0).
           cimg_forXY(img,x,y) {
               FMState state = grid[img.width()*(img.height()-y-1)+x].getState();
               unsigned int val = 0;
               if (state == FMState::NARROW)
                   val = 127;
               else if (state == FMState::OPEN)
                   val = 255;
               img(x,y) = val;
           }
           //img.map( CImg<double>::jet_LUT256() );
           name += "FMStates";
           img.display(name.c_str(), false);
       }

};

#endif /* GRIDPLOTTER_H_ */
