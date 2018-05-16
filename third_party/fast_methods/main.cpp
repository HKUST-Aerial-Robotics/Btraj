/* n-dimensional Fast Marching example with the main functions used */

#include <iostream>
#include <cmath>
#include <chrono>
#include <array>
#include <string>

#include "fmm/fmdata/fmcell.h"
#include "ndgridmap/ndgridmap.hpp"
#include "console/console.h"
#include "fmm/fmm.hpp"
#include "fm2/fm2.hpp"
#include "fm2/fm2star.hpp"
#include "fmm/fmdata/fmfibheap.hpp"
#include "fmm/fmdata/fmpriorityqueue.hpp"
#include "fmm/fmdata/fmdaryheap.hpp"
#include "fmm/fmdata/fmdaryheap.hpp"
#include "io/maploader.hpp"
#include "io/gridplotter.hpp"
#include "io/gridwriter.hpp"
#include "io/gridpoints.hpp"
#include "gradientdescent/gradientdescent.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char ** argv)
{
    constexpr unsigned int ndims = 2; // Setting two dimensions.
    constexpr unsigned int ndims3 = 3; // Setting three dimensions.

    time_point<std::chrono::steady_clock> start, end; // Time measuring.
    double time_elapsed;

    console::info("Parsing input arguments.");
    string filename1, filename2, filename_vels;
    console::parseArguments(argc,argv, "-map1", filename1);
    console::parseArguments(argc,argv, "-map2", filename2);
    console::parseArguments(argc,argv, "-vel", filename_vels);

    console::info("Creating grid from image and testing Fast Marching Method..");
    nDGridMap<FMCell, ndims> grid;
    MapLoader::loadMapFromImg(filename2.c_str(), grid);

    std::array<unsigned int, ndims> coords_init, coords_goal;
    GridPoints::selectMapPoints(grid, coords_init, coords_goal);

    vector<unsigned int> init_points;
    unsigned int idx, goal;
    grid.coord2idx(coords_init, idx);
    init_points.push_back(idx);
    grid.coord2idx(coords_goal, goal);

    FMM< nDGridMap<FMCell, ndims> > fmm;
    fmm.setEnvironment(&grid);
    fmm.setInitialAndGoalPoints(init_points, goal);
    fmm.compute();
        cout << "\tElapsed FM time: " << fmm.getTime() << " ms" << endl;
    console::info("Plotting the results and saving into test_fm.txt");
    GridPlotter::plotArrivalTimes(grid);
    GridWriter::saveGridValues("test_fm.txt", grid);

    console::info("Computing gradient descent ");
    typedef typename std::vector< std::array<double, ndims> > Path; // A bit of short-hand.

    Path path;
    std::vector <double> path_velocity; // Velocities profile

        start = steady_clock::now();
    GradientDescent< nDGridMap<FMCell, ndims> > grad;
    grad.apply(grid,goal,path,path_velocity);
        end = steady_clock::now();
        time_elapsed = duration_cast<milliseconds>(end-start).count();
        cout << "\tElapsed gradient descent time: " << time_elapsed << " ms" << endl;
    GridWriter::savePath("test_path.txt", grid, path);
    GridWriter::savePathVelocity("path_velocity.txt", grid, path, path_velocity);
    GridPlotter::plotMapPath(grid,path);

    console::info("Testing Fast Marching Square Method.");
    path_velocity.clear();
    grid.coord2idx(coords_goal, goal);
    Path pathFM2;

    MapLoader::loadMapFromImg(filename2.c_str(), grid);

    FM2<nDGridMap<FMCell, ndims>> fm2;
    fm2.setEnvironment(&grid);
    fm2.setInitialAndGoalPoints(init_points, goal);
    fm2.compute();
        cout << "\tElapsed FM2 time: " << fm2.getTime() << " ms" << endl;
        start = steady_clock::now();
    fm2.computePath(&pathFM2, &path_velocity);
        end = steady_clock::now();
        time_elapsed = duration_cast<milliseconds>(end-start).count();
        cout << "\tElapsed gradient descent time: " << time_elapsed << " ms" << endl;

    GridPlotter::plotOccupancyPath(grid, pathFM2, fm2.getName());
    GridPlotter::plotMapPath(grid, pathFM2, fm2.getName());

    console::info("Testing FM2*.");
    Path pathFM2star;
    FM2Star<nDGridMap<FMCell, ndims>> fm2star;
    fm2star.setEnvironment(&grid);
    fm2star.setInitialAndGoalPoints(init_points, goal);
    fm2star.compute();
        cout << "\tElapsed FM2* time: " << fm2star.getTime() << " ms" << endl;
        start = steady_clock::now();
    fm2star.computePath(&pathFM2star, &path_velocity);
        end = steady_clock::now();
        time_elapsed = duration_cast<milliseconds>(end-start).count();
        cout << "\tElapsed gradient descent time: " << time_elapsed << " ms" << endl;

    console::info("Comparing FM2, FM2* paths.");
    std::vector <Path> paths;
    paths.push_back(pathFM2);
    paths.push_back(pathFM2star);

    GridPlotter::plotMapPaths(grid,paths);

    console::info("Now let's try different velocities.");
    nDGridMap<FMCell, ndims> grid_vels;
    MapLoader::loadMapFromImg(filename_vels.c_str(), grid_vels);

    GridPlotter::plotOccupancyMap(grid_vels);

    FMM< nDGridMap<FMCell, ndims> , FMFibHeap<> > fmm_vels;
    init_points.clear();
    init_points.push_back(80000); // Just an init point as any other.
    fmm_vels.setEnvironment(&grid_vels);
    fmm_vels.setInitialPoints(init_points);
    fmm_vels.compute();
        cout << "\tElapsed FM time: " << fmm_vels.getTime() << " ms" << endl;

    console::info("Plotting the results ");
    GridPlotter::plotArrivalTimes(grid_vels);

    console::info("Testing 3D!");
    nDGridMap<FMCell, ndims3> grid3 (std::array<unsigned int,ndims3>{100,100,50});
    init_points.clear();
    grid3.coord2idx(std::array<unsigned int,ndims3>{50,50,25},idx); // Reusing the previous int.
    init_points.push_back(idx);
    grid3.coord2idx(std::array<unsigned int, ndims3> {20, 10, 45}, goal);

    FMM< nDGridMap<FMCell, ndims3> > fmm3;
    fmm3.setEnvironment(&grid3);
    fmm3.setInitialAndGoalPoints(init_points, goal);
    fmm3.compute();
        cout << "\tElapsed FM time: " << fmm3.getTime() << " ms" << endl;

    console::info("Saving into file test_fm3d.txt");
    GridWriter::saveGridValues("test_fm3d.txt", grid3);

    console::info("Testing 3D gradient descent.");
    typedef typename std::vector< std::array<double, ndims3> > Path3D; // A bit of short-hand.

    grid3.coord2idx(std::array<unsigned int, ndims3> {20, 10, 45}, goal);

    Path3D path3D;
        start = steady_clock::now();
    GradientDescent< nDGridMap<FMCell, ndims3> > grad3D;
    grad3D.apply(grid3,goal,path3D,path_velocity);
        end = steady_clock::now();
        time_elapsed = duration_cast<milliseconds>(end-start).count();
        cout << "\tElapsed gradient descent time: " << time_elapsed << " ms" << endl;
    GridWriter::savePath("test_path3d.txt", grid3, path3D);

    return 0;
}
