/*! \class BenchmarkCFG
    \brief Configures a Benchmark class from a CFG file. 
    Inspired in the OMPLapp benchmarking software: http://ompl.kavrakilab.org

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

#ifndef BENCHMARKCFG_HPP_
#define BENCHMARKCFG_HPP_

#include <string>
#include <unordered_map>
#include <algorithm>

#include <boost/program_options.hpp>
#include <boost/algorithm/string/case_conv.hpp>

#include "benchmark.hpp"
#include "../io/maploader.hpp"

#include "../fmm/fmm.hpp"
#include "../fmm/sfmm.hpp"
#include "../fmm/fmmstar.hpp"
#include "../fmm/sfmmstar.hpp"
#include "../fmm/fim.hpp"
#include "../fmm/gmm.hpp"
#include "../fmm/ufmm.hpp"
#include "../fmm/fsm.hpp"
#include "../fmm/lsm.hpp"
#include "../fmm/ddqm.hpp"

/// \todo the getter functions do not check if the types are admissible.
/// \todo does not have support for multiple starts or goals.
/// \todo the way ctor parameters are given could be improved (as declareParams in OMPL - command pattern).

class BenchmarkCFG {

    public:

        /** \brief Requires readOptions to be manually called after this constructor. */
        BenchmarkCFG
        () {}

        /** \brief Creates an object from a CFG file. Automatically calls readOptions. */
        BenchmarkCFG
        (const char * filename) {
            readOptions(filename);
        }

        /** \brief Parses the CFG file given. */
        bool readOptions(const char * filename)
        {
            static const std::vector<std::string> knownSolvers = {
                "fmm", "fmmstar", "fmmfib", "fmmfibstar", "sfmm", "sfmmstar",
                "gmm", "fim", "ufmm", "fsm", "lsm", "ddqm" // Add solver here.
            };

            std::fstream cfg(filename);
            if (!cfg.good())
            {
               std::string error("Unable to open file: ");
               error += filename;
               console::error(error);
               return 0;
            }

            path_ = boost::filesystem::absolute( boost::filesystem::path(filename) );
            boost::filesystem::path name = path_.filename();
            name.replace_extension("");

            boost::program_options::options_description desc;
            desc.add_options()
                ("grid.file",          boost::program_options::value<std::string>(),                             "Path to load a velocities map from image.")
                ("grid.text",          boost::program_options::value<std::string>(),                             "Path to load a velocities map from a .grid file.")
                ("grid.ndims",         boost::program_options::value<std::string>()->default_value("2"),         "Number of dimensions.")
                ("grid.cell",          boost::program_options::value<std::string>()->default_value("FMCell"),    "Type of cell. FMCell by default.")
                ("grid.dimsize",       boost::program_options::value<std::string>()->default_value("200,200"),   "Size of dimensions: N,M,O...")
                ("grid.leafsize",      boost::program_options::value<std::string>()->default_value("1"),         "Leafsize (assuming cubic cells).")
                ("problem.start",      boost::program_options::value<std::string>()->required(),                 "Start point: s1,s2,s3...")
                ("problem.goal",       boost::program_options::value<std::string>()->default_value("nan"),       "Goal point: g1,g2,g3... By default no goal point.")
                ("benchmark.name",     boost::program_options::value<std::string>()->default_value(name.string()), "Name of the benchmark.")
                ("benchmark.runs",     boost::program_options::value<std::string>()->default_value("10"),        "Number of runs per solver.")
                ("benchmark.savegrid", boost::program_options::value<std::string>()->default_value("0"),         "Save grid values of each run.");

            boost::program_options::variables_map vm;
            boost::program_options::parsed_options po = boost::program_options::parse_config_file(cfg, desc, true);
            boost::program_options::store(po, vm);
            boost::program_options::notify(vm);
            cfg.close();

            // Storing registered options, intended for grid, problem and benchmark configuration.
            for (const auto& var : vm)
                options_[var.first] = boost::any_cast<std::string>(var.second.value());

            // Analyze unregistered options, intended for solvers.
            std::vector<std::string> unr = boost::program_options::collect_unrecognized(po.options, boost::program_options::exclude_positional);
            for (std::size_t i = 0 ; i < unr.size()/2 ; ++i)
            {
                std::string key = boost::to_lower_copy(unr[i*2]);
                std::string val = unr[i*2 + 1];

                if (key.substr(0,8) == "solvers.")
                {
                    std::string solver = key.substr(8);
                    // If it is a known solver, save it together with its constructor
                    // parameters (if given).
                    for (std::size_t i = 0; i < knownSolvers.size(); ++i)
                        if (solver == knownSolvers[i]) {
                            solverNames_.push_back(knownSolvers[i]);
                            ctorParams_.push_back(val);
                        }
                }
            }

            return 1;
       }

        template <class grid_t, class cell_t>
        void configure
        (Benchmark<grid_t> & b) {
            b.setSaveLog(true);
            b.setName(getValue<std::string>("benchmark.name"));
            b.setSaveGrid(getValue<unsigned int>("benchmark.savegrid"));
            b.setNRuns(getValue<unsigned int>("benchmark.runs"));
            b.setPath(boost::filesystem::path("results"));
            b.fromCFG(true);

            for (size_t i = 0; i < solverNames_.size(); ++i)
            {
                const std::string & name = solverNames_[i];

                bool defaultCtor = true;
                if (!ctorParams_[i].empty())
                    defaultCtor = false;

                Solver<grid_t> * solver;

                if (defaultCtor) {
                    if (name == "fmm")
                        solver = new FMM<grid_t>();
                    else if (name == "fmmstar")
                        solver = new FMMStar<grid_t>();
                    else if (name == "fmmfib")
                        solver = new FMM<grid_t, FMFibHeap<cell_t> >("FMMFib");
                    else if (name == "fmmfibstar")
                        solver = new FMMStar<grid_t,  FMFibHeap<cell_t> >("FMMFib*");
                    else if (name == "sfmm")
                        solver = new SFMM<grid_t>("SFMM");
                    else if (name == "sfmmstar")
                        solver = new SFMMStar<grid_t>("SFMM*");
                    else if (name == "gmm")
                        solver = new GMM<grid_t>();
                    else if (name == "fim")
                        solver = new FIM<grid_t>();
                    else if (name == "ufmm")
                        solver = new UFMM<grid_t>();
                    else if (name == "fsm")
                        solver = new FSM<grid_t>();
                    else if (name == "lsm")
                        solver = new LSM<grid_t>();
                    else if (name == "ddqm")
                        solver = new DDQM<grid_t>();
                    // Add solver here.

                    else
                        continue;
                }
                else { // Create solvers with specified constructor parameters.
                    std::vector<std::string> p(split(ctorParams_[i]));

                    // FMM and FMM*
                    if (name == "fmm")
                        solver = new FMM<grid_t>(ctorParams_[i].c_str());
                    else if (name == "fmmstar") {
                        if (p.size() == 1)
                            solver = new FMMStar<grid_t>(p[0].c_str());
                        else if (p.size() == 2) {
                            if (p[1] == "TIME")
                                solver = new FMMStar<grid_t>(p[0].c_str(), TIME);
                            else if (p[1] == "DISTANCE")
                                solver = new FMMStar<grid_t>(p[0].c_str(), DISTANCE);
                        }
                    }
                    // FMMFib and FMMFib*
                    else if (name == "fmmfib")
                        solver = new FMM<grid_t, FMFibHeap<cell_t> >(ctorParams_[i].c_str());
                    else if (name == "fmmfibstar") {
                        if (p.size() == 1)
                            solver = new FMMStar<grid_t, FMFibHeap<cell_t>>(p[0].c_str());
                        else if (p.size() == 2) {
                            if (p[1] == "TIME")
                                solver = new FMMStar<grid_t, FMFibHeap<cell_t>>(p[0].c_str(), TIME);
                            else if (p[1] == "DISTANCE")
                                solver = new FMMStar<grid_t, FMFibHeap<cell_t>>(p[0].c_str(), DISTANCE);
                        }
                    }
                    // SFMM and SFMM*
                    else if (name == "sfmm")
                        solver = new SFMM<grid_t, cell_t>(ctorParams_[i].c_str());
                    else if (name == "sfmmstar") {
                        if (p.size() == 1)
                            solver = new SFMMStar<grid_t, cell_t>(p[0].c_str());
                        else if (p.size() == 2) {
                            if (p[1] == "TIME")
                                solver = new SFMMStar<grid_t, cell_t>(p[0].c_str(), TIME);
                            else if (p[1] == "DISTANCE")
                                solver = new SFMMStar<grid_t, cell_t>(p[0].c_str(), DISTANCE);
                        }
                    }
                    // GMM
                    else if (name == "gmm") {
                        if (p.size() == 1)
                            solver = new GMM<grid_t>(p[0].c_str());
                        else if (p.size() == 2)
                            solver = new GMM<grid_t>(p[0].c_str(), boost::lexical_cast<double>(p[1]));
                    }
                    // FIM
                    else if (name == "fim") {
                        if (p.size() == 1)
                            solver = new FIM<grid_t>(p[0].c_str());
                        else if (p.size() == 2)
                            solver = new FIM<grid_t>(p[0].c_str(), boost::lexical_cast<double>(p[1]));
                    }
                    // UFMM
                    else if (name == "ufmm") {
                        if (p.size() == 1)
                            solver = new UFMM<grid_t>(p[0].c_str());
                        else if (p.size() == 2)
                            solver = new UFMM<grid_t>(p[0].c_str(), boost::lexical_cast<unsigned>(p[1]));
                        else if (p.size() == 3)
                            solver = new UFMM<grid_t>(p[0].c_str(), boost::lexical_cast<unsigned>(p[1]), boost::lexical_cast<double>(p[2]));
                    }
                    // FSM
                    else if (name == "fsm") {
                        if (p.size() == 1)
                            solver = new FSM<grid_t>(p[0].c_str());
                        else if (p.size() == 2)
                            solver = new FSM<grid_t>(p[0].c_str(), boost::lexical_cast<unsigned>(p[1]));
                    }
                    // LSM
                    else if (name == "lsm") {
                        if (p.size() == 1)
                            solver = new LSM<grid_t>(p[0].c_str());
                        else if (p.size() == 2)
                            solver = new LSM<grid_t>(p[0].c_str(), boost::lexical_cast<unsigned>(p[1]));
                    }
                    // DDQM
                    else if (name == "ddqm") {
                        solver = new LSM<grid_t>(p[0].c_str());
                    }
                    // Add solver here.

                    else
                        continue;
                }

                b.addSolver(solver);
            }

            constexpr size_t N = grid_t::getNDims();
            grid_t * grid = new grid_t();
            if (options_.find("grid.file") != options_.end())
                MapLoader::loadMapFromImg(options_.find("grid.file")->second.c_str(), *grid);
            else if (options_.find("grid.text") != options_.end()) {
                if(!MapLoader::loadMapFromText(options_.find("grid.text")->second.c_str(), *grid))
                    exit(1);
            }
            else {
                const std::string & strToSplit = options_.find("grid.dimsize")->second;
                std::array<unsigned int, N> dimSize = splitAndCast<unsigned int, N>(strToSplit);
                grid->resize(dimSize);
            }

            grid->setLeafSize(getValue<double>("grid.leafsize"));

            const std::string & strToSplit2 = options_.find("problem.start")->second;
            std::array<unsigned int, N> startCoords = splitAndCast<unsigned int, N>(strToSplit2);
            unsigned int startIdx;
            std::vector<unsigned int> startIndices;
            grid->coord2idx(startCoords, startIdx);
            startIndices.push_back(startIdx);

            const std::string & strToSplit3 = options_.find("problem.goal")->second;
            if (strToSplit3 != "nan")
            {
                unsigned int goalIdx;
                std::array<unsigned int, N> goalCoords = splitAndCast<unsigned int, N>(strToSplit3);
                grid->coord2idx(goalCoords, goalIdx);
                b.setInitialAndGoalPoints(startIndices, goalIdx);
            }
            else
                b.setInitialPoints(startIndices);

            b.setEnvironment(grid);
        }

        /** \brief Get the value for a given key (option). */
        template<typename T>
        T getValue
        (const std::string & key) const {
            std::unordered_map<std::string, std::string>::const_iterator iter = options_.find(key);
            if (iter != options_.end())
                return boost::lexical_cast<T>(iter->second);
            else
                return boost::lexical_cast<T>(0);
        }

    private:
        // Based on http://stackoverflow.com/a/236803/2283531
        /** \brief From a string of format XXX,YY,ZZZ,... splis the N elements and cast as type T (comma-separated) as an array. */
        template <typename T, size_t N>
        std::array<T,N> splitAndCast
        (const std::string & s)
        {
            std::array<T,N> elems;
            std::stringstream ss(s);
            std::string item;
            for (size_t i = 0; i < N; ++i)
            {
                std::getline(ss, item, ',');
                elems[i] = boost::lexical_cast<T>(item);
            }
            return elems;
        }

        /** \brief From a string of format XXX,YY,ZZZ,... splis the elements (comma-separated) as a vector of strings. */
        std::vector<std::string> split
        (const std::string & s) {
            std::vector<std::string> elems;
            std::stringstream ss(s);
            std::string item;
            size_t pos = 1, str_size = s.length()-1;
            while (pos < str_size) {
                std::getline(ss, item, ',');
                elems.push_back(item);
                pos += item.length();
            }
            return elems;
        }

        /** \brief Stores the names of the parsed solvers. */
        std::vector<std::string> solverNames_;
        
        /** \brief Stores the constructor parameters for the parsed solvers. */
        std::vector<std::string> ctorParams_;
        
        /** \brief Option-value map.*/
        std::unordered_map<std::string, std::string> options_;

        /** \brief Stores the absolute path to the CFG file. */
        boost::filesystem::path path_;
};

#endif /* BENCHMARKCFG_HPP_*/
