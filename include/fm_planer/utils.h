#ifndef _UTILS_
#define _UTILS_

#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/SDF.h>
#include <sdf_tools/sdf.hpp>
#include <sdf_tools/collision_map.hpp>
#include <sdf_tools/dynamic_spatial_hashed_collision_map.hpp>

#include "../third_party/fast_marching/fmm/fmdata/fmcell.h"
#include "../third_party/fast_marching/ndgridmap/ndgridmap.hpp"
#include "../third_party/fast_marching/fmm/fmm.hpp"
#include "../third_party/fast_marching/fmm/ufmm.hpp"
#include "../third_party/fast_marching/fmm/gmm.hpp"
#include "../third_party/fast_marching/fmm/fmmstar.hpp"
#include "../third_party/fast_marching/fmm/sfmm.hpp"
#include "../third_party/fast_marching/fmm/sfmmstar.hpp"
#include "../third_party/fast_marching/gradientdescent/gradientdescent.hpp"
//#include <array>

#define _PI M_PI
typedef nDGridMap<FMCell, 3> FMGrid3D;
typedef array<unsigned int, 3> Coord3D;
typedef typename std::vector< array<double, 3> > Path3D; 

#endif