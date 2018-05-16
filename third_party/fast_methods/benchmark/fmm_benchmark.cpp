/*! \brief Automatically configures and runs a Benchmark from a CFG file.

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

#include <iostream>
#include <cmath>
#include <array>
#include <string>

#include <boost/variant.hpp>

#include "../fmm/fmdata/fmcell.h"
#include "../ndgridmap/ndgridmap.hpp"

#include "../fmm/fmm.hpp"
#include "../fmm/fim.hpp"
#include "../fmm/gmm.hpp"
#include "../fmm/ufmm.hpp"

#include "../fmm/fmdata/fmfibheap.hpp"
#include "../fmm/fmdata/fmdaryheap.hpp"
#include "../fmm/fmdata/fmpriorityqueue.hpp"

#include "../benchmark/benchmark.hpp"
#include "../benchmark/benchmarkcfg.hpp"

using namespace std;

int main(int argc, const char ** argv)
{
    // Parse input.
    if (argc < 2)
    {
        std::cerr << "Usage:\n\t " << argv[0] << " problem.cfg" << std::endl;
        return 1;
    }

    // Parse the CFG file.
    BenchmarkCFG bcfg;
    if (bcfg.readOptions(argv[1]))
    {
        // If FMCell is used...
        if(bcfg.getValue<std::string>("grid.cell") == "FMCell")
        {
            // ... and for dimensions...
            switch (bcfg.getValue<unsigned int>("grid.ndims"))
            {
                case 2:
                {
                    Benchmark<nDGridMap<FMCell,2> > b;
                    bcfg.configure<nDGridMap<FMCell,2>, FMCell>(b);
                    b.run();
                    break;
                }
                case 3:
                {
                    Benchmark<nDGridMap<FMCell,3> > b;
                    bcfg.configure<nDGridMap<FMCell,3>, FMCell>(b);
                    b.run();
                    break;
                }
                // Include here new dimensions copying, pasting and changing x.
                /*case x:
                {
                    Benchmark<nDGridMap<FMCell,x> > b;
                    bcfg.configure<nDGridMap<FMCell,x>, FMCell>(b);
                    b.run();
                    break;
                }*/
            }
        }
        else // else if (bcfg.getValue<std::string>("grid.cell") == "MyCell") 
        {
            // Include here new celltypes and include the corresponding switch dimensions as for FMCell.
            // ... and for dimensions...
            /*switch (bcfg.getValue<unsigned int>("grid.ndims"))
            {
                case 2:
                {
                    Benchmark<nDGridMap<MyCell,2> > b;
                    bcfg.configure<nDGridMap<MyCell,2>, MyCell>(b);
                    b.run();
                    break;
                }
                //...
            }*/
        }
    }
}
