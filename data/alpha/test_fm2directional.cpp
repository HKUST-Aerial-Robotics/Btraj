/* n dimensional Fast Marching example with the main functions used */

#include <iostream>
#include <cmath>
#include <chrono>
#include <array>
#include <string>
#include <algorithm>

#include "../fmm/fmdata/fmdirectionalcell.h"
#include "../ndgridmap/ndgridmap.hpp"
#include "../console/console.h"
#include "../fm2/fm2dir.hpp"
#include "../fmm/fmdata/fmfibheap.hpp"
#include "../fmm/fmdata/fmpriorityqueue.hpp"
#include "../fmm/fmdata/fmdaryheap.hpp"
#include "../fmm/fmdata/fmdaryheap.hpp"
#include "../io/maploader.hpp"

using namespace std;
using namespace std::chrono;

int main(int argc, const char ** argv)
{
    constexpr int ndims2 = 2; // Setting two dimensions.

    console::info("Parsing input arguments.");
    string filename;
    if (argc > 2)
        console::parseArguments(argc,argv, "-map", filename);
    else {
        console::info("No enough arguments given. Loading default example map: examples/data/grid.txt");
        filename = "../data/grid.txt";
    }

    // A bit of shorthand.
    typedef nDGridMap<FMDirectionalCell, ndims2> FMGrid2D;
    typedef array<int, ndims2> Coord2D;

    time_point<std::chrono::system_clock> start, end; // Time measuring.
    double time_elapsed;

    FMGrid2D grid_fm2;

    if(!MapLoaderText::loadMapFromText(filename.c_str(), grid_fm2))
        exit(1);

    Coord2D init_point = {377, 664};
    Coord2D goal_point = {379, 91};
    vector<int> init_points;
    int idx, goal;
    grid_fm2.coord2idx(init_point , idx);
    init_points.push_back(idx);
    grid_fm2.coord2idx(goal_point , goal);

    std::vector<Solver<FMGrid2D>*> solvers;
    solvers.push_back(new FM2Dir<FMGrid2D>("FM2Dir_Dary"));
    solvers.push_back(new FM2Dir<FMGrid2D, FMFibHeap<FMCell> >("FM2Dir_Fib"));
    solvers.push_back(new FM2Dir<FMGrid2D, FMPriorityQueue<FMCell> >("FM2Dir_SFMM"));

    for (Solver<FMGrid2D>* s :solvers)
    {
        s->setEnvironment(&grid_fm2);
            start = system_clock::now();
        s->setInitialAndGoalPoints(init_points, goal_idx);
        s->compute();
            end = system_clock::now();
            time_elapsed = duration_cast<milliseconds>(end-start).count();
            cout << "\tElapsed "<< s->getName() <<" time: " << time_elapsed << " ms" << '\n';
        GridPlotter::plotArrivalTimes(grid_fm2);
    }

    // Preventing memory leaks.
    for (auto & s : solvers)
        delete s;
}
