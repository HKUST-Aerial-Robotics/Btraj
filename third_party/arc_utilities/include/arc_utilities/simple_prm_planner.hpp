#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <functional>
#include <chrono>
#include <random>
#include <omp.h>
#ifdef ENABLE_PARALLEL_ROADMAP
    #ifndef ENABLE_PARALLEL_K_NEAREST_NEIGHBORS
        #define ENABLE_PARALLEL_K_NEAREST_NEIGHBORS
        #include <arc_utilities/arc_helpers.hpp>
        #undef ENABLE_PARALLEL_K_NEAREST_NEIGHBORS
    #else
        #include <arc_utilities/arc_helpers.hpp>
    #endif
#else
    #include <arc_utilities/arc_helpers.hpp>
#endif
#include <arc_utilities/dijkstras.hpp>


#ifndef SIMPLE_PRM_PLANNER_HPP
#define SIMPLE_PRM_PLANNER_HPP

namespace simple_prm_planner
{
    class SimpleGeometricPrmPlanner
    {
    protected:

        SimpleGeometricPrmPlanner() {}

    public:

        template<typename T, typename Allocator=std::allocator<T>>
        static arc_dijkstras::Graph<T, Allocator> BuildRoadMap(const std::function<T(void)>& sampling_fn,
                                                               const std::function<double(const T&, const T&)>& distance_fn,
                                                               const std::function<bool(const T&)>& state_validity_check_fn,
                                                               const std::function<bool(const T&, const T&)>& edge_validity_check_fn,
                                                               const std::function<bool(void)>& termination_check_fn,
                                                               const size_t K, const bool distance_is_symmetric=true)
        {
            std::function<double(const arc_dijkstras::GraphNode<T, Allocator>&, const T&)> graph_distance_fn = [&] (const arc_dijkstras::GraphNode<T, Allocator>& node, const T& state) { return distance_fn(node.GetValueImmutable(), state); };
            arc_dijkstras::Graph<T, Allocator> roadmap;
            while (!termination_check_fn())
            {
                const T random_state = sampling_fn();
                if (!state_validity_check_fn(random_state))
                {
                    continue;
                }
                const std::vector<std::pair<int64_t, double>> nearest_neighbors = arc_helpers::GetKNearestNeighbors(roadmap.GetNodesImmutable(), random_state, graph_distance_fn, K);
                const int64_t new_node_index = roadmap.AddNode(random_state);
                // Parallelize the collision-checking and distance computation
                std::vector<std::pair<double, double>> nearest_neighbors_distances(nearest_neighbors.size());
#ifdef ENABLE_PARALLEL_ROADMAP
                #pragma omp parallel for
#endif
                for (size_t idx = 0; idx < nearest_neighbors.size(); idx++)
                {
                    const std::pair<int64_t, double>& nearest_neighbor = nearest_neighbors[idx];
                    const int64_t nearest_neighbor_index = nearest_neighbor.first;
                    const double nearest_neighbor_distance = nearest_neighbor.second;
                    const T& nearest_neighbor_state = roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable();
                    if (edge_validity_check_fn(nearest_neighbor_state, random_state))
                    {
                        if (distance_is_symmetric)
                        {
                            nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, nearest_neighbor_distance);
                        }
                        else
                        {
                            const double reverse_distance = distance_fn(random_state, roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable());
                            nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, reverse_distance);
                        }
                    }
                    else
                    {
                        nearest_neighbors_distances[idx] = std::make_pair(-1.0, -1.0);
                    }
                }
                // THIS MUST BE SERIAL - add edges to roadmap
                for (size_t idx = 0; idx < nearest_neighbors.size(); idx++)
                {
                    const std::pair<int64_t, double>& nearest_neighbor = nearest_neighbors[idx];
                    const int64_t nearest_neighbor_index = nearest_neighbor.first;
                    const std::pair<double, double>& nearest_neighbor_distances = nearest_neighbors_distances[idx];
                    if (nearest_neighbor_distances.first >= 0.0 && nearest_neighbor_distances.second >= 0.0)
                    {
                        roadmap.AddEdgeBetweenNodes(nearest_neighbor_index, new_node_index, nearest_neighbor_distances.first);
                        roadmap.AddEdgeBetweenNodes(new_node_index, nearest_neighbor_index, nearest_neighbor_distances.second);
                    }
                }
            }
            return roadmap;
        }

        template<typename T, typename Allocator=std::allocator<T>>
        static void UpdateRoadMapEdges(arc_dijkstras::Graph<T, Allocator>& roadmap, const std::function<bool(const T&, const T&)>& edge_validity_check_fn, const std::function<double(const T&, const T&)>& distance_fn)
        {
            assert(roadmap.CheckGraphLinkage());
#ifdef ENABLE_PARALLEL_ROADMAP
            #pragma omp parallel for
#endif
            for (size_t current_node_index = 0; current_node_index < roadmap.GetNodesImmutable().size(); current_node_index++)
            {
                arc_dijkstras::GraphNode<T, Allocator>& current_node = roadmap.GetNodeMutable(current_node_index);
                std::vector<arc_dijkstras::GraphEdge>& current_node_out_edges = current_node.GetOutEdgesMutable();
                for (size_t out_edge_idx = 0; out_edge_idx < current_node_out_edges.size(); out_edge_idx++)
                {
                    arc_dijkstras::GraphEdge& current_out_edge = current_node_out_edges[out_edge_idx];
                    const int64_t other_node_idx = current_out_edge.GetToIndex();
                    arc_dijkstras::GraphNode<T, Allocator>& other_node = roadmap.GetNodeMutable(other_node_idx);
                    std::vector<arc_dijkstras::GraphEdge>& other_node_in_edges = other_node.GetInEdgesMutable();
                    const double updated_weight = (edge_validity_check_fn(current_node.GetValueImmutable(), other_node.GetValueImmutable())) ? distance_fn(current_node.GetValueImmutable(), other_node.GetValueImmutable()) : std::numeric_limits<double>::infinity();
                    // Update our out edge
                    current_out_edge.SetWeight(updated_weight);
                    // Update the other node's in edges
                    for (size_t in_edge_idx = 0; in_edge_idx < other_node_in_edges.size(); in_edge_idx++)
                    {
                        arc_dijkstras::GraphEdge& other_in_edge = other_node_in_edges[in_edge_idx];
                        if (other_in_edge.GetFromIndex() == current_node_index)
                        {
                            other_in_edge.SetWeight(updated_weight);
                        }
                    }
                }
            }
        }

        template<typename T, typename Allocator=std::allocator<T>>
        static std::pair<std::vector<T, Allocator>, double> QueryPathAndAddNodes(const T& start, const T& goal, arc_dijkstras::Graph<T, Allocator>& roadmap, const std::function<bool(const T&, const T&)>& edge_validity_check_fn, const std::function<double(const T&, const T&)>& distance_fn, const size_t K, const bool distance_is_symmetric=true)
        {
            // Add the start node to the roadmap
            std::function<double(const arc_dijkstras::GraphNode<T, Allocator>&, const T&)> start_distance_fn = [&] (const arc_dijkstras::GraphNode<T, Allocator>& node, const T& state) { return distance_fn(state, node.GetValueImmutable()); };
            const std::vector<std::pair<int64_t, double>> start_nearest_neighbors = arc_helpers::GetKNearestNeighbors(roadmap.GetNodesImmutable(), start, start_distance_fn, K);
            const int64_t start_node_index = roadmap.AddNode(start);
            // Parallelize the collision-checking and distance computation
            std::vector<std::pair<double, double>> start_nearest_neighbors_distances(start_nearest_neighbors.size());
#ifdef ENABLE_PARALLEL_ROADMAP
            #pragma omp parallel for
#endif
            for (size_t idx = 0; idx < start_nearest_neighbors.size(); idx++)
            {
                const std::pair<int64_t, double>& nearest_neighbor = start_nearest_neighbors[idx];
                const int64_t nearest_neighbor_index = nearest_neighbor.first;
                const double nearest_neighbor_distance = nearest_neighbor.second;
                const T& nearest_neighbor_state = roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable();
                if (edge_validity_check_fn(nearest_neighbor_state, start))
                {
                    if (distance_is_symmetric)
                    {
                        start_nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, nearest_neighbor_distance);
                    }
                    else
                    {
                        const double reverse_distance = distance_fn(roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable(), start);
                        start_nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, reverse_distance);
                    }
                }
                else
                {
                    start_nearest_neighbors_distances[idx] = std::make_pair(-1.0, -1.0);
                }
            }
            // THIS MUST BE SERIAL - add edges to roadmap
            for (size_t idx = 0; idx < start_nearest_neighbors.size(); idx++)
            {
                const std::pair<int64_t, double>& nearest_neighbor = start_nearest_neighbors[idx];
                const int64_t nearest_neighbor_index = nearest_neighbor.first;
                const std::pair<double, double>& nearest_neighbor_distances = start_nearest_neighbors_distances[idx];
                if (nearest_neighbor_distances.first >= 0.0 && nearest_neighbor_distances.second >= 0.0)
                {
                    roadmap.AddEdgeBetweenNodes(start_node_index, nearest_neighbor_index, nearest_neighbor_distances.first);
                    roadmap.AddEdgeBetweenNodes(nearest_neighbor_index, start_node_index, nearest_neighbor_distances.second);
                }
            }
            // Add the goal node to the roadmap
            std::function<double(const arc_dijkstras::GraphNode<T, Allocator>&, const T&)> goal_distance_fn = [&] (const arc_dijkstras::GraphNode<T, Allocator>& node, const T& state) { return distance_fn(node.GetValueImmutable(), state); };
            const std::vector<std::pair<int64_t, double>> goal_nearest_neighbors = arc_helpers::GetKNearestNeighbors(roadmap.GetNodesImmutable(), goal, goal_distance_fn, K);
            const int64_t goal_node_index = roadmap.AddNode(goal);
            // Parallelize the collision-checking and distance computation
            std::vector<std::pair<double, double>> goal_nearest_neighbors_distances(goal_nearest_neighbors.size());
#ifdef ENABLE_PARALLEL_ROADMAP
            #pragma omp parallel for
#endif
            for (size_t idx = 0; idx < goal_nearest_neighbors.size(); idx++)
            {
                const std::pair<int64_t, double>& nearest_neighbor = goal_nearest_neighbors[idx];
                const int64_t nearest_neighbor_index = nearest_neighbor.first;
                const double nearest_neighbor_distance = nearest_neighbor.second;
                const T& nearest_neighbor_state = roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable();
                if (edge_validity_check_fn(nearest_neighbor_state, goal))
                {
                    if (distance_is_symmetric)
                    {
                        goal_nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, nearest_neighbor_distance);
                    }
                    else
                    {
                        const double reverse_distance = distance_fn(goal, roadmap.GetNodeImmutable(nearest_neighbor_index).GetValueImmutable());
                        goal_nearest_neighbors_distances[idx] = std::make_pair(nearest_neighbor_distance, reverse_distance);
                    }
                }
                else
                {
                    goal_nearest_neighbors_distances[idx] = std::make_pair(-1.0, -1.0);
                }
            }
            // THIS MUST BE SERIAL - add edges to roadmap
            for (size_t idx = 0; idx < goal_nearest_neighbors.size(); idx++)
            {
                const std::pair<int64_t, double>& nearest_neighbor = goal_nearest_neighbors[idx];
                const int64_t nearest_neighbor_index = nearest_neighbor.first;
                const std::pair<double, double>& nearest_neighbor_distances = goal_nearest_neighbors_distances[idx];
                if (nearest_neighbor_distances.first >= 0.0 && nearest_neighbor_distances.second >= 0.0)
                {
                    roadmap.AddEdgeBetweenNodes(nearest_neighbor_index, goal_node_index, nearest_neighbor_distances.first);
                    roadmap.AddEdgeBetweenNodes(goal_node_index, nearest_neighbor_index, nearest_neighbor_distances.second);
                }
            }
            // Call Dijkstra's
            const auto dijkstras_solution = arc_dijkstras::SimpleDijkstrasAlgorithm<T, Allocator>::PerformDijkstrasAlgorithm(roadmap, goal_node_index);
            // Extract solution
            const std::pair<std::vector<int64_t>, std::vector<double>>& solution_map_distances = dijkstras_solution.second;
            const double start_node_distance = solution_map_distances.second[start_node_index];
            if (std::isinf(start_node_distance))
            {
                return std::make_pair(std::vector<T, Allocator>(), std::numeric_limits<double>::infinity());
            }
            else
            {
                std::vector<int64_t> solution_path_indices;
                solution_path_indices.push_back(start_node_index);
                int64_t previous_index = solution_map_distances.first[start_node_index];
                while (previous_index >= 0)
                {
                    const int64_t current_index = previous_index;
                    solution_path_indices.push_back(current_index);
                    if (current_index == goal_node_index)
                    {
                        break;
                    }
                    else
                    {
                        previous_index = solution_map_distances.first[current_index];
                    }
                }
                std::vector<T, Allocator> solution_path;
                solution_path.reserve(solution_path_indices.size());
                for (size_t idx = 0; idx < solution_path_indices.size(); idx++)
                {
                    const int64_t solution_path_index = solution_path_indices[idx];
                    const T& solution_path_state = roadmap.GetNodeImmutable(solution_path_index).GetValueImmutable();
                    solution_path.push_back(solution_path_state);
                }
                solution_path.shrink_to_fit();
                return std::make_pair(solution_path, start_node_distance);
            }
        }

        template<typename T, typename Allocator=std::allocator<T>>
        static std::pair<std::vector<T, Allocator>, double> QueryPath(const T& start, const T& goal, const arc_dijkstras::Graph<T, Allocator>& roadmap, const std::function<bool(const T&, const T&)>& edge_validity_check_fn, const std::function<double(const T&, const T&)>& distance_fn, const size_t K, const bool distance_is_symmetric=true)
        {
            arc_dijkstras::Graph<T, Allocator> working_copy = roadmap;
            return QueryPathAndAddNodes(start, goal, working_copy, edge_validity_check_fn, distance_fn, K, distance_is_symmetric);
        }
    };
}

#endif // SIMPLE_PRM_PLANNER
