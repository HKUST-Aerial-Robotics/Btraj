/*! \class GradientDescent
    \brief Implements gradient descent to be used together with nDGridMap
    and FMCell or derived types.
   
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

#ifndef GRADIENTDESCENT_H_
#define GRADIENTDESCENT_H_

#include <iostream>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <numeric>

#include "../ndgridmap/ndgridmap.hpp"
#include "../fm/fmdata/fmcell.h"

template <typename T>
T sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <class grid_t> class GradientDescent {

    /** \brief Shorthand for number of dimensions. */
    static constexpr size_t ndims_ = grid_t::getNDims();

    /** \brief Shorthand for coordinates. */
    typedef typename std::array<unsigned int, ndims_> Coord;

    /** \brief Shorhand for real points. */
    typedef typename std::array<double, ndims_> Point;

    /** \brief Shorthand for path type of real points. */
    typedef typename std::vector <Point> Path;

    public:
        /** \brief Computes the path over grid from idx to a minimum and extracts the velocity in every point.

           Simple gradient approximation is used: dimension 0: gx = f((x-1,y)+f(x+1,y))/2
           dimension 1: gy = f((x,y-1)+f(x,y+1))/2
           and so on. */

        /*Note: This function is the original function that retrives the path.
                It has some problems: when using FM with heuristic, often results in a bad path (zig-zag), or can not reach the target.
        */

        static void apply
        (grid_t & grid, unsigned int & idx, Path & path, std::vector <double> & path_velocity, std::vector <double> & time, double step = 1) 
        { 
          int iter_num = 0;
          Coord current_coord;
          Point current_point;
          Coord dimsize = grid.getDimSizes();

          std::array<unsigned int, ndims_-1> d_; //  Same as nDGridMap class auxiliar array d_.
          //cout<<"check d_[]"<<endl;
          d_[0] = dimsize[0];
          //cout<<d_[0]<<endl;
          for (size_t i = 1; i < ndims_; ++i)
          {
              d_[i] = dimsize[i]*d_[i-1];
              //cout<<d_[i]<<endl;
          }

          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to int.
          path.push_back(current_point);
          path_velocity.push_back(grid[idx].getVelocity());
          time.push_back(grid[idx].getArrivalTime());
          
          std::array<double, ndims_> grads;
          //std::array<double, ndims_> grads_tmp;

          while(grid[idx].getArrivalTime() != 0) 
          {
              if(iter_num >= 30000)
              {
                cout<<"can not find a path by gradient descent"<<endl;
                break;
              }

              iter_num ++;
              // Every iteration the gradient is computed for all dimensions. If is infinite, we convert it to 1 (keeping the sign).
              // The static_cast are necessary because the conversion between coordinate (we check the value in coordinates) and points
              // (the path is composed by continuous points).

              grads[0] =   - grid[idx-1].getValue()/2 + grid[idx+1].getValue()/2;
                          
              if (std::isinf(grads[0]) )
                  grads[0] = 0.5 * sgn<double>(grads[0]);

              if (std::isnan(grads[0]))
                  grads[0] = 0.0;

              size_t idx_1 = min(max(idx-d_[0], (unsigned int)0), d_[2] - 1);
              size_t idx_2 = min(max(idx+d_[0], (unsigned int)0), d_[2] - 1);

              grads[1] = - grid[idx_1].getValue()/2 + grid[idx_2].getValue()/2;
                          
              if (std::isinf(grads[1]))
                  grads[1] = 0.5 * sgn<double>(grads[1]);
              if (std::isnan(grads[1]))
                  grads[1] = 0.0;

              idx_1 = min(max(idx-d_[1], (unsigned int)0), d_[2] - 1);
              idx_2 = min(max(idx+d_[1], (unsigned int)0), d_[2] - 1);
              grads[2] = - grid[idx_1].getValue()/2 + grid[idx_2].getValue()/2;

              if (std::isinf(grads[2]))
                  grads[2] = 0.5 * sgn<double>(grads[2]);

              if (std::isnan(grads[2]))
                  grads[2] = 0.0;
              double grad_norm = sqrt(grads[0] * grads[0] + grads[1] * grads[1] + grads[2] * grads[2]);

              // Updating points
              for (size_t i = 0; i < ndims_; ++i) 
              {
                  // Moving the point in dim i.
                  current_point[i] = current_point[i] - step*grads[i]/grad_norm;
                  current_coord[i] = int(current_point[i] +0.5);
              }
              path.push_back(current_point);
              path_velocity.push_back(grid[idx].getVelocity());
              grid.coord2idx(current_coord,idx);
              time.push_back(grid[idx].getArrivalTime());
          }
          //Adding exactly the last point at the end.
          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to double.
          path.push_back(current_point);
          path_velocity.push_back(grid[idx].getVelocity());
          time.push_back(grid[idx].getArrivalTime());
      }

      /*Note: This function is implemented by myself, the most differences are:
              1. use 9 neighbours of each voxel to calculate the gradient in one direction (except for neighbor with a inf value)
              2. if the obtained gradient is ill (inf, NaN or all 0), locally search a neighbour with biggest value drop to go.
      */

      static int extract_path
      (grid_t & grid, unsigned int & idx, Path & path, std::vector <double> & path_velocity, std::vector <double> & time, double step = 1) 
      { 
          int iter_num = 0;
          Coord current_coord;
          Point current_point;
          Coord dimsize = grid.getDimSizes();

          std::array<unsigned int, ndims_-1> d_; //  Same as nDGridMap class auxiliar array d_.
          d_[0] = dimsize[0];
          
          for (size_t i = 1; i < ndims_; ++i)
              d_[i] = dimsize[i]*d_[i-1];

          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to int.
          path.push_back(current_point);
          time.push_back(grid[idx].getArrivalTime());
          
          path_velocity.push_back(grid[idx].getVelocity());
          
          std::array<double, ndims_> grads;
          unsigned int index1, index2;

          while(grid[idx].getArrivalTime() != 0) 
          {
              if(iter_num >= 10000)
              {
                cout<<"can not find a path by gradient descent"<<endl;
                return -1;
              }
              iter_num ++;
              
              // get X dimension gradient
              // ++++++++++++++++++++++++++++++++++++++++++
              //cout<<"start an iteration"<<endl;
              {
                  int idx_x_bias_l = (int)idx - (int)1;
                  int idx_x_bias_u = (int)idx + (int)1;

                  index1   = min( (unsigned int)max(idx_x_bias_l, 0), d_[2] - 2);
                  index2   = min( (unsigned int)max(idx_x_bias_u, 0), d_[2] - 2);

                  double grad = 0.0;
                  int cnt = 9;
                  // get surrounding 8 points
                  for(int i = -1; i < 2; i++)
                  {
                    for(int j = -1; j < 2; j++)
                    { 
                        int index_n_l = (int)index1 + i * (int)d_[0] + j * (int)d_[1];
                        unsigned int index_s_l = min( (unsigned int)max( index_n_l, 0), d_[2] - 2);
                        double value_l = grid[index_s_l].getValue();

                        int index_n_u = (int)index2 + i * (int)d_[0] + j * (int)d_[1];
                        unsigned int index_s_u = min( (unsigned int)max(index_n_u, 0), d_[2] - 2);
                        double value_u = grid[index_s_u].getValue();
                        
                        if( std::isinf(value_l) || std::isinf(value_u))
                            cnt--;
                        else
                            grad += (value_u - value_l) / 2.0;
                    }
                  }
                  grad /= cnt;
                  grads[0] = grad;
              }
              // get Y dimension gradient
              // **************************************************************** //    
              {
                  int idx_y_bias_l = (int)idx - (int)d_[0];
                  int idx_y_bias_u = (int)idx + (int)d_[0];

                  index1 = min( (unsigned int)max(idx_y_bias_l, 0), d_[2] - 2);
                  index2 = min( (unsigned int)max(idx_y_bias_u, 0), d_[2] - 2);

                  double grad = 0.0;
                  int cnt = 9;
                  // get surrounding 8 points
              
                  for(int i = -1; i < 2; i++)
                  {
                    for(int j = -1; j < 2; j++)
                    { 
                        int index_n_l = (int)index1 + i * (int)1 + j * (int)d_[1];
                        unsigned int index_s_l = min( (unsigned int)max( index_n_l, 0), d_[2] - 2);
                        double value_l = grid[index_s_l].getValue();

                        int index_n_u = (int)index2 + i * (int)1 + j * (int)d_[1];
                        unsigned int index_s_u = min( (unsigned int)max( index_n_u, 0), d_[2] - 2);
                        double value_u = grid[index_s_u].getValue();

                        if( std::isinf(value_l) || std::isinf(value_u))
                            cnt--;
                        else
                            grad += (value_u - value_l) / 2.0;
                    }
                  }
                  grad /= cnt;
                  grads[1] = grad;
              }
              
              // get Z dimension gradient
              // **************************************************************** //    

              //cout<<"get z grad"<<endl;
              {
                  int idx_z_bias_l = (int)idx - (int)d_[1];
                  int idx_z_bias_u = (int)idx + (int)d_[1];

                  index1 = min( (unsigned int)max(idx_z_bias_l, 0), d_[2] - 2);
                  index2 = min( (unsigned int)max(idx_z_bias_u, 0), d_[2] - 2);

                  double grad = 0.0;
                  int cnt = 9;
                  // get surrounding 8 points
                  for(int i = -1; i < 2; i++)
                  {
                    for(int j = -1; j < 2; j++)
                    {   
                        int index_n_l = (int)index1 + i * (int)1 + j * (int)d_[0];
                        unsigned int index_s_l = min( (unsigned int)max( index_n_l, 0), d_[2] - 2);
                        double value_l = grid[index_s_l].getValue();

                        int index_n_u = (int)index2 + i * (int)1 + j * (int)d_[0];
                        unsigned int index_s_u = min( (unsigned int)max( index_n_u, 0), d_[2] - 2);
                        double value_u = grid[index_s_u].getValue();

                        if( std::isinf(value_l) || std::isinf(value_u))
                            cnt--;
                        else
                            grad += (value_u - value_l) / 2.0;
                    }
                  }
                  grad /= cnt;
                  grads[2] = grad;
              }

              double max_dif = -1000.0;
              double dif;
              unsigned int best_idx = idx;

              if( std::isinf(grads[0]) || std::isinf(grads[1]) || std::isinf(grads[2])
                  || std::isnan(grads[0]) || std::isnan(grads[1]) || std::isnan(grads[2])
                  || (grads[0] == 0.0 && grads[1] == 0.0 && grads[2] == 0.0) )
              {   
                  for(int i = -1; i < 2; i++)
                    for(int j = -1; j < 2; j++)
                      for(int k = -1; k < 2; k++)
                      {
                          if( i == 0 && j == 0 && k == 0)
                              continue;

                          int idx_delta = i + (int)d_[0] * j + (int)d_[1] * k;
                          int idx_n = (int)idx + idx_delta;

                          unsigned int idx_tmp = (unsigned int)idx_n;
                          idx_tmp = min(max(idx_tmp, (unsigned int)0), d_[2] - 2);

                          if( std::isinf(grid[idx_tmp].getValue()))
                              continue;

                          dif = (-grid[idx_tmp].getValue() + grid[idx].getValue()) / sqrt( double(i*i + j*j + k*k) );
                          if(dif > max_dif)
                          {
                              max_dif = dif;
                              best_idx = idx_tmp;
                          }
                      }

                  best_idx = min(max(best_idx, (unsigned int)0), d_[2] - 2);
                  grid.idx2coord(best_idx, current_coord);
                  std::copy_n( current_coord.begin(), ndims_, current_point.begin() );
              } 
              else
              {   
                  double grad_norm = sqrt(grads[0] * grads[0] + grads[1] * grads[1] + grads[2] * grads[2]);
                  // Updating points
                  for (size_t i = 0; i < ndims_; ++i) 
                  {
                      current_point[i] = current_point[i] - step*grads[i]/grad_norm; 
                      current_coord[i] = int(current_point[i] + 0.5);
                  }
              }

              path.push_back(current_point);
              path_velocity.push_back(grid[idx].getVelocity());
              time.push_back(grid[idx].getArrivalTime());

              grid.coord2idx(current_coord,idx);
          }

          //Adding exactly the last point at the end.
          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to double.
          path.push_back(current_point);
          time.push_back(grid[idx].getArrivalTime());
          path_velocity.push_back(grid[idx].getVelocity());

          return 1;
      }

      static int gradient_descent
      (grid_t & grid, unsigned int & idx, Path & path, std::vector <double> & path_velocity, std::vector <double> & time, double step = 1) 
      { 
          int iter_num = 0;
          Coord current_coord;
          Point current_point;
          Coord dimsize = grid.getDimSizes();

          std::array<unsigned int, ndims_-1> d_; //  Same as nDGridMap class auxiliar array d_.
          d_[0] = dimsize[0];
          
          for (size_t i = 1; i < ndims_; ++i)
              d_[i] = dimsize[i]*d_[i-1];

          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to int.
          path.push_back(current_point);
          time.push_back(grid[idx].getArrivalTime());
          
          path_velocity.push_back(grid[idx].getVelocity());
          
          while(grid[idx].getArrivalTime() != 0) 
          {
              if(iter_num >= 10000)
              {
                cout<<"can not find a path by gradient descent"<<endl;
                return -1;
              }
              iter_num ++;

              double max_dif = -10000.0;
              double dif;
              unsigned int best_idx = idx;

              for(int i = -1; i < 2; i++)
                for(int j = -1; j < 2; j++)
                  for(int k = -1; k < 2; k++)
                  {
                      if( i == 0 && j == 0 && k == 0)
                          continue;

                      int idx_delta = i + (int)d_[0] * j + (int)d_[1] * k;
                      int idx_n = (int)idx + idx_delta;

                      unsigned int idx_tmp = (unsigned int)idx_n;
                      idx_tmp = min(max(idx_tmp, (unsigned int)0), d_[2] - 2);

                      if( std::isinf(grid[idx_tmp].getValue()))
                          continue;

                      dif = (-grid[idx_tmp].getValue() + grid[idx].getValue()) / sqrt( double(i*i + j*j + k*k) );
                      if(dif > max_dif)
                      {
                          max_dif = dif;
                          best_idx = idx_tmp;
                      }
                  }

              best_idx = min(max(best_idx, (unsigned int)0), d_[2] - 2);
              grid.idx2coord(best_idx, current_coord);
              std::copy_n( current_coord.begin(), ndims_, current_point.begin() );

              path.push_back(current_point);
              path_velocity.push_back(grid[idx].getVelocity());
              time.push_back(grid[idx].getArrivalTime());

              grid.coord2idx(current_coord,idx);
          }

          //Adding exactly the last point at the end.
          grid.idx2coord(idx, current_coord);
          std::copy_n( current_coord.begin(), ndims_, current_point.begin() ); // Cast to double.
          path.push_back(current_point);
          time.push_back(grid[idx].getArrivalTime());
          path_velocity.push_back(grid[idx].getVelocity());

          return 1;
      }
};

#endif /* GRADIENTDESCENT_H_*/