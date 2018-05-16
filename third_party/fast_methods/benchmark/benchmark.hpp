/*! \class Benchmark
    \brief This class provides the utilities to benchmark Fast Marching Solvers (configuration, running and logging).
    It works for FMM (any heap and SFMM), FIM and UFMM. It has not been tested
    with FM2 solvers.
    
    By default, it will save a log file in a generated folder called results.
    
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

#ifndef BENCHMARK_HPP_
#define BENCHMARK_HPP_

#include <chrono>
#include <limits>

#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include "../fmm/solver.hpp"
#include "../io/gridwriter.hpp"

template <class grid_t>
class Benchmark {

    public:

        Benchmark
        (unsigned int saveGrid = 0, bool saveLog = true) :
        saveGrid_(saveGrid),
        saveLog_(saveLog),
        runID_(0),
        nruns_(10),
        path_("results"),
        name_("benchmark"),
        fromCFG_(false){}

        virtual ~Benchmark()
        {
            clear();
        }

        /** \brief Set true if the benchmark is configured from a cfg file. */
        void fromCFG
        (bool cfg) {
            fromCFG_ = cfg;
        }

        /** \brief Adds a new solver to be benchmarked. */
        void addSolver
        (Solver<grid_t>* solver) {
            solvers_.push_back(solver);
        }

        /** \brief Sets the saveGrid_ flag. */
        void setSaveGrid
        (unsigned int s) {
            saveGrid_ = s;
        }

        /** \brief  Sets the environment (grid map) the benchmark will be run on. */
        void setEnvironment
        (grid_t* grid) {
            grid_ = grid;
        }

        /** \brief  Set number of runs for each solver. */
        void setNRuns
        (unsigned int n) {
            nruns_ = n;
        }

        /** \brief  Set the path where results will be saved. */
        void setPath
        (const boost::filesystem::path & path) {
            path_ = path;
        }

        /** \brief Sets the flag to save the log to a file (true) or output in terminal (false). */
        void setSaveLog
        (bool s) {
            saveLog_ = s;
        }

        /** \brief Sets the initial and goal points (indices) for the solvers. */
        void setInitialAndGoalPoints
        (const std::vector<unsigned int> & init_points, unsigned int goal_idx) {
            init_points_ = init_points;
            goal_idx_ = goal_idx;
        }

        /** \brief Sets the initial points (indices for the solvers). */
        void setInitialPoints(const std::vector<unsigned int> & init_points)
        {
            setInitialAndGoalPoints(init_points, -1);
        }

        /** \brief Automatically runs all the solvers. */
        void run
        () {
            if (saveGrid_)
            {
                boost::filesystem::create_directory(path_);
                boost::filesystem::create_directory(path_ / name_);
            }
            else if (saveLog_)
                boost::filesystem::create_directory(path_);

            boost::progress_display showProgress (solvers_.size()*nruns_);

            configSolvers();

            logConfig();

            for (Solver<grid_t>* s :solvers_)
            {
                for (unsigned int i = 0; i < nruns_; ++i)
                {
                    ++runID_;
                    s->reset();
                    s->compute();
                    logRun(s);

                    if (saveGrid_ == 2)
                        saveGrid(s);

                    ++showProgress;
                }
                if (saveGrid_ == 1)
                    saveGrid(s);
                s->reset();                
            }

            if (saveLog_)
                saveLog();
            else {
                console::info("Benchmark log format:");
                std::cout << "Name\t#Runs\t#Dims\tDim1...DimN\t#Starts\tStartIdx\tGoalIdx"<<'\n';
                std::cout << "RunID\tName\tTime (ms)" << '\n';
                std::cout << log_.str() << '\n';
            }
        }

        /** \brief  Logs the last run of solver s. */
        void logRun
        (const Solver<grid_t>* s)
        {
            std::ios init(NULL);
            init.copyfmt(std::cout);
            formatID();
            log_ << '\n' << fmtID_;

            std::cout.copyfmt(init);
            log_ << '\t' << s->getName() << "\t" << s->getTime();
        }

        /** \brief Saves the grid values result of the last run of solver s. */
        void saveGrid
        (Solver<grid_t>* s) const {
            thread_local boost::filesystem::path filename;
            if (saveGrid_ == 1)
                filename = path_ / name_ / s->getName();
            if (saveGrid_ == 2)
                filename = path_ / name_ / fmtID_;
            
            filename.replace_extension(".grid");
            GridWriter::saveGridValues(filename.string().c_str(), *(s->getGrid()));
        }

        /** \brief Saves the log to a file: benchmark_name.log */
        void saveLog
        () const {
            std::ofstream ofs (path_.string() + "/" + name_ + ".log");
            ofs << log_.rdbuf();
            ofs.close();
        }

        /** \brief Sets the name of the benchmark. */
        void setName
        (const std::string & n)
        {
            name_ = n;
        }

        /** \brief Clears and free memory allocated by the benchmark. */
        void clear
        () {
            for (auto & s : solvers_)
                delete s;

            solvers_.clear();

            if(fromCFG_)
                delete grid_;
        }

    private:
        /** \brief Logs the benchmark configuration. */
        void logConfig
        () {
            log_ << name_ << '\t' << nruns_ << '\t' << grid_->getNDims()
                 << '\t' << grid_->getDimSizesStr();

            // Saving starting points.
            log_ << init_points_.size() << '\t';
            std::copy(init_points_.begin(), init_points_.end(), std::ostream_iterator<unsigned int>(log_, "\t"));

            if (int(goal_idx_) == -1)
                log_ << "nan";
            else
                log_ << goal_idx_;
        }

        /** \brief Configures each solver so they are ready to be run. */
        void configSolvers
        () {
            for (Solver<grid_t>* s :solvers_)
            {
                s->setEnvironment(grid_);
                s->setInitialAndGoalPoints(init_points_, goal_idx_);
            }
        }

        /** \brief Formats as a string the run ID. */
        void formatID
        () {
            thread_local std::ostringstream oss; // To optimize a bit.
            oss.str("");
            oss.clear();
            oss << std::setw(4) << std::setfill('0') << runID_;
            fmtID_ = oss.str();
        }

        /** \brief Stores the solvers to be run in the benchmark. */
        std::vector<Solver<grid_t>*>                        solvers_;
        
        /** \brief Grid the benchmark will be run on. */
        grid_t *                                            grid_;

        /** \brief Indices of the initial points. */
        std::vector<unsigned int>                           init_points_;
        
        /** \brief Index of the goal point. */
        unsigned int                                        goal_idx_;

        /** \brief Time measurement variables. */        
        std::chrono::time_point<std::chrono::system_clock>  start_, end_;

        /** \brief Log stream. */
        std::stringstream                                   log_;

        /** \brief If 1, the resulting grids (times) of the first run of each solver are saved to files.
             If 2, the grids for all runs of each solver are saved. */
        unsigned int                                        saveGrid_;
        
        /** \brief  If true, the log is saved to file. Output on terminal otherwise. */
        bool                                                saveLog_;

        /** \brief ID of the current run. */
        unsigned int                                        runID_;
        
        /** \brief Formatted run ID: */
        std::string                                         fmtID_;
        
        /** \brief Number of runs for each solver. */
        unsigned int                                        nruns_;

        /** \brief Path were the results will be saved. */
        boost::filesystem::path                             path_;
        
        /** \brief Benchmark name. */
        std::string                                         name_;

        /** \brief If true, benchmark configured from CFG file, used to selectively free memory. */
        bool                                                fromCFG_;
};

#endif /* BENCHMARK_HPP_*/
