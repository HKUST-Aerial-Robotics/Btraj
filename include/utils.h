#ifndef _UTILS_
#define _UTILS_

#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/SDF.h>
#include <sdf_tools/sdf.hpp>
#include <sdf_tools/collision_map.hpp>
#include <sdf_tools/dynamic_spatial_hashed_collision_map.hpp>

#include "../third_party/fast_methods/gradientdescent/gradientdescent.hpp"
#include "../third_party/fast_methods/ndgridmap/ndgridmap.hpp"
#include "../third_party/fast_methods/fm/fmdata/fmcell.h"
#include "../third_party/fast_methods/fm/fmm.hpp"
#include "../third_party/fast_methods/fm/fmmstar.hpp"

/*
#include "../third_party/fast_methods/fm/ufmm.hpp"
#include "../third_party/fast_methods/fm/gmm.hpp"
#include "../third_party/fast_methods/fm/fim.hpp"
#include "../third_party/fast_methods/fm/fsm.hpp"
#include "../third_party/fast_methods/fm/lsm.hpp"
#include "../third_party/fast_methods/fm/ddqm.hpp"
#include "../third_party/fast_methods/fm/sfmm.hpp"
#include "../third_party/fast_methods/fm/sfmmstar.hpp"
*/

#define _PI M_PI
typedef nDGridMap<FMCell, 3> FMGrid3D;
typedef array<unsigned int, 3> Coord3D;
typedef typename std::vector< array<double, 3> > Path3D; 

#endif