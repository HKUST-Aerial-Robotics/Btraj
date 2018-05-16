#ifndef FIRST_ORDER_VISIBILITY_GRAPH_HPP
#define FIRST_ORDER_VISIBILITY_GRAPH_HPP

#include <memory>
#include <functional>
#include <Eigen/Core>
#include <omp.h>
#include <ros/ros.h>
#include <visualization_msgs/MarkerArray.h>
#include "arc_utilities/arc_helpers.hpp"

namespace arc_utilities
{
    class FirstOrderVisibilityGraph
    {
        public:
            typedef std::pair<ssize_t, ssize_t> ConfigType;
            typedef std::pair<ConfigType, double> ConfigAndDistType;
            typedef std::function<bool(const ssize_t row, const ssize_t col)> ValidityCheckFnType;

            static double ConfigTypeDistance(const ConfigType& c1, const ConfigType& c2)
            {
                return Eigen::Vector2d((double)(c1.first - c2.first), (double)(c1.second - c2.second)).norm();
            }

            struct BestFirstSearchComparator
            {
                public:
                    // Defines a "less" operation"; by using "greater" then the smallest element will appear at the top of the priority queue
                    bool operator()(const ConfigAndDistType& c1, const ConfigAndDistType& c2) const
                    {
                        // If expected distances are different, we want to explore the one with the smaller expected distance
                        return (c1.second > c2.second);
                    }
            };

            static bool CheckFirstOrderVisibility(const ssize_t rows, const ssize_t cols, const ValidityCheckFnType& validity_check_fn, const bool visualization_enabled = true)
            {
                assert(rows > 0 && cols > 0);
                typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> ArrayXb;

                ros::NodeHandle nh;
                ros::Publisher marker_pub = nh.advertise<visualization_msgs::Marker>("visualization_marker", 1000, true);

                visualization_msgs::Marker marker;
                {
                    marker.header.frame_id = "mocap_world";
                    marker.type = visualization_msgs::Marker::POINTS;
                    marker.action = visualization_msgs::Marker::ADD;
                    marker.ns = "explored_states";
                    marker.id = 1;
                    marker.scale.x = 1.0;
                    marker.scale.y = 1.0;
                    marker.header.stamp = ros::Time::now();
                    marker_pub.publish(marker);
                    marker.color = arc_helpers::RGBAColorBuilder<std_msgs::ColorRGBA>::MakeFromFloatColors(0.0, 0.0, 1.0, 1.0);
                }

                const ConfigType start(0, 0), goal(rows - 1, cols - 1);
                const auto heuristic_distance_fn = [&goal] (const ConfigType& config) { return ConfigTypeDistance(config, goal); };

                std::priority_queue<ConfigAndDistType, std::vector<ConfigAndDistType>, BestFirstSearchComparator> frontier;
                ArrayXb explored = ArrayXb::Constant(rows, cols, false);

                frontier.push(ConfigAndDistType(start, heuristic_distance_fn(start)));

                bool path_found = false;
                std::cout << "Entering explore loop\n\n";
                while (!path_found && ros::ok() && frontier.size() > 0)
                {
                    const ConfigAndDistType current = frontier.top();
                    frontier.pop();
                    const ConfigType& current_config = current.first;

                    // Visualization code
                    if (visualization_enabled)
                    {
                        geometry_msgs::Point p;
                        p.x = current_config.first;
                        p.y = current_config.second;
                        p.z = 0;
                        ++marker.id;
                        marker.points.push_back(p);

                        if (marker.id % 1000 == 0)
                        {
                            marker.header.stamp = ros::Time::now();
                            marker_pub.publish(marker);
                            marker.points.clear();
                            usleep(10);
                        }
                    }

                    if (current_config.first == goal.first && current_config.second == goal.second)
                    {
                        if (visualization_enabled)
                        {
                            std::cout << "Reached goal!\n";
                            std::cout << PrettyPrint::PrettyPrint(current_config, true, " ") << std::endl << std::flush;
                            std::cout << std::endl;
                            marker.header.stamp = ros::Time::now();
                            marker_pub.publish(marker);
                        }
                        path_found = true;
                    }
                    // Double check if we've already explored this node:
                    //    a single node can be inserted into the frontier multiple times at the same or different priorities
                    //    so we want to avoid the expense of re-exploring it, and just discard this one once we pop it
                    else if (explored(current_config.first, current_config.second) == false)
                    {
                        explored(current_config.first, current_config.second) = true;

                        // Expand the node to find all neighbours, adding them to the frontier if we have not already explored them
                        const auto neighbours = GetNeighbours(current_config, rows, cols, validity_check_fn);
                        for (const auto neighbour : neighbours)
                        {
                            // Check if we've already explored this neighbour to avoid re-adding it to the frontier
                            if (explored(neighbour.first, neighbour.second) == false)
                            {
                                frontier.push(ConfigAndDistType(neighbour, heuristic_distance_fn(neighbour)));
                            }
                        }
                    }
                }

                if (visualization_enabled)
                {
                    marker.header.stamp = ros::Time::now();
                    marker_pub.publish(marker);
                }

                return path_found;
            }

        private:
            FirstOrderVisibilityGraph() {}

            static std::vector<ConfigType> GetNeighbours(const ConfigType& config, const ssize_t rows, const ssize_t cols, const ValidityCheckFnType& validity_check_fn)
            {
                std::vector<ConfigType> neighbours;
                neighbours.reserve(8);

                const ssize_t row_min = std::max(0L, config.first - 1);
                const ssize_t row_max = std::min(rows - 1, config.first + 1);

                const ssize_t col_min = std::max(0L, config.second - 1);
                const ssize_t col_max = std::min(cols - 1, config.second + 1);

                for (ssize_t col = col_min; col <= col_max; col++)
                {
                    for (ssize_t row = row_min; row <= row_max; row++)
                    {
                        if (!(row == config. first && col == config.second) && validity_check_fn(row, col) == true)
                        {
                            neighbours.push_back(ConfigType(row, col));
                        }
                    }
                }

                return neighbours;
            }

    };
}

#endif // FIRST_ORDER_VISIBILITY_GRAPH_HPP
