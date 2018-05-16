#ifndef DIJKSTRAS_HPP
#define DIJKSTRAS_HPP

#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <Eigen/Geometry>
#include <arc_utilities/arc_helpers.hpp>

namespace arc_dijkstras
{
    class GraphEdge
    {
        protected:

            int64_t from_index_;
            int64_t to_index_;
            double weight_;

        public:

            static uint64_t Serialize(const GraphEdge& edge, std::vector<uint8_t>& buffer)
            {
                return edge.SerializeSelf(buffer);
            }

            static std::pair<GraphEdge, uint64_t> Deserialize(const std::vector<uint8_t>& buffer, const uint64_t current)
            {
                GraphEdge temp_edge;
                const uint64_t bytes_read = temp_edge.DeserializeSelf(buffer, current);
                return std::make_pair(temp_edge, bytes_read);
            }

            GraphEdge(const int64_t from_index, const int64_t to_index, const double weight)
                : from_index_(from_index), to_index_(to_index), weight_(weight)
            {}

            GraphEdge()
                : from_index_(-1), to_index_(-1), weight_(0.0)
            {}

            uint64_t SerializeSelf(std::vector<uint8_t>& buffer) const
            {
                const uint64_t start_buffer_size = buffer.size();
                arc_helpers::SerializeFixedSizePOD<int64_t>(from_index_, buffer);
                arc_helpers::SerializeFixedSizePOD<int64_t>(to_index_, buffer);
                arc_helpers::SerializeFixedSizePOD<double>(weight_, buffer);
                // Figure out how many bytes were written
                const uint64_t end_buffer_size = buffer.size();
                const uint64_t bytes_written = end_buffer_size - start_buffer_size;
                return bytes_written;
            }

            uint64_t DeserializeSelf(const std::vector<uint8_t>& buffer, const uint64_t current)
            {
                assert(current < buffer.size());
                uint64_t current_position = current;
                const std::pair<int64_t, uint64_t> deserialized_from_index = arc_helpers::DeserializeFixedSizePOD<int64_t>(buffer, current_position);
                from_index_ = deserialized_from_index.first;
                current_position += deserialized_from_index.second;
                const std::pair<int64_t, uint64_t> deserialized_to_index = arc_helpers::DeserializeFixedSizePOD<int64_t>(buffer, current_position);
                to_index_ = deserialized_to_index.first;
                current_position += deserialized_to_index.second;
                const std::pair<double, uint64_t> deserialized_weight = arc_helpers::DeserializeFixedSizePOD<double>(buffer, current_position);
                weight_ = deserialized_weight.first;
                current_position += deserialized_weight.second;
                // Figure out how many bytes were read
                const uint64_t bytes_read = current_position - current;
                return bytes_read;
            }

            bool operator==(const GraphEdge& other) const
            {
                return (from_index_ == other.GetFromIndex() && to_index_ == other.GetToIndex() && weight_ == other.GetWeight());
            }

            std::string Print() const
            {
                return std::string("(" + std::to_string(from_index_) + "->" + std::to_string(to_index_) + ") : " + std::to_string(weight_));
            }

            int64_t GetFromIndex() const
            {
                return from_index_;
            }

            int64_t GetToIndex() const
            {
                return to_index_;
            }

            double GetWeight() const
            {
                return weight_;
            }

            void SetFromIndex(const int64_t new_from_index)
            {
                from_index_ = new_from_index;
            }

            void SetToIndex(const int64_t new_to_index)
            {
                to_index_ = new_to_index;
            }

            void SetWeight(const double new_weight)
            {
                weight_ = new_weight;
            }
    };

    inline std::ostream& operator<< (std::ostream& stream, const GraphEdge& edge)
    {
        stream << edge.Print();
        return stream;
    }

    template<typename NodeValueType, typename Allocator=std::allocator<NodeValueType>>
    class GraphNode
    {
        protected:

            NodeValueType value_;
            double distance_;
            std::vector<GraphEdge> in_edges_;
            std::vector<GraphEdge> out_edges_;

        public:

            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            static uint64_t Serialize(const GraphNode<NodeValueType, Allocator>& node, std::vector<uint8_t>& buffer, const std::function<uint64_t(const NodeValueType&, std::vector<uint8_t>&)>& value_serializer)
            {
                return node.SerializeSelf(buffer, value_serializer);
            }

            static std::pair<GraphNode<NodeValueType, Allocator>, uint64_t> Deserialize(const std::vector<uint8_t>& buffer, const uint64_t current, const std::function<std::pair<NodeValueType, uint64_t>(const std::vector<uint8_t>&, const uint64_t)>& value_deserializer)
            {
                GraphNode<NodeValueType, Allocator> temp_node;
                const uint64_t bytes_read = temp_node.DeserializeSelf(buffer, current, value_deserializer);
                return std::make_pair(temp_node, bytes_read);
            }

            GraphNode(const NodeValueType& value, const double distance, const std::vector<GraphEdge>& new_in_edges, const std::vector<GraphEdge>& new_out_edges)
                : value_(value), distance_(distance), in_edges_(new_in_edges), out_edges_(new_out_edges)
            {}

            GraphNode(const NodeValueType& value)
                : value_(value), distance_(std::numeric_limits<double>::infinity())
            {}

            GraphNode()
                : distance_(std::numeric_limits<double>::infinity())
            {}

            uint64_t SerializeSelf(std::vector<uint8_t>& buffer, const std::function<uint64_t(const NodeValueType&, std::vector<uint8_t>&)>& value_serializer) const
            {
                const uint64_t start_buffer_size = buffer.size();
                // Serialize the value
                value_serializer(value_, buffer);
                // Serialize the distance
                arc_helpers::SerializeFixedSizePOD<double>(distance_, buffer);
                // Serialize the in edges
                arc_helpers::SerializeVector<GraphEdge>(in_edges_, buffer, GraphEdge::Serialize);
                // Serialize the in edges
                arc_helpers::SerializeVector<GraphEdge>(out_edges_, buffer, GraphEdge::Serialize);
                // Figure out how many bytes were written
                const uint64_t end_buffer_size = buffer.size();
                const uint64_t bytes_written = end_buffer_size - start_buffer_size;
                return bytes_written;
            }

            uint64_t DeserializeSelf(const std::vector<uint8_t>& buffer, const uint64_t current, const std::function<std::pair<NodeValueType, uint64_t>(const std::vector<uint8_t>&, const uint64_t)>& value_deserializer)
            {
                assert(current < buffer.size());
                uint64_t current_position = current;
                // Deserialize the value
                const std::pair<NodeValueType, uint64_t> value_deserialized = value_deserializer(buffer, current_position);
                value_ = value_deserialized.first;
                current_position += value_deserialized.second;
                // Deserialize the distace
                const std::pair<double, uint64_t> distance_deserialized = arc_helpers::DeserializeFixedSizePOD<double>(buffer, current_position);
                distance_ = distance_deserialized.first;
                current_position += distance_deserialized.second;
                // Deserialize the in edges
                const std::pair<std::vector<GraphEdge>, uint64_t> in_edges_deserialized = arc_helpers::DeserializeVector<GraphEdge>(buffer, current_position, GraphEdge::Deserialize);
                in_edges_ = in_edges_deserialized.first;
                current_position += in_edges_deserialized.second;
                // Deserialize the out edges
                const std::pair<std::vector<GraphEdge>, uint64_t> out_edges_deserialized = arc_helpers::DeserializeVector<GraphEdge>(buffer, current_position, GraphEdge::Deserialize);
                out_edges_ = out_edges_deserialized.first;
                current_position += out_edges_deserialized.second;
                // Figure out how many bytes were read
                const uint64_t bytes_read = current_position - current;
                return bytes_read;
            }

            std::string Print() const
            {
                std::ostringstream strm;
                strm << "Node : " << distance_ << " In Edges : ";
                if (in_edges_.size() > 0)
                {
                    strm << in_edges_[0].Print();
                    for (size_t idx = 1; idx < in_edges_.size(); idx++)
                    {
                        strm << ", " << in_edges_[idx].Print();
                    }
                }
                strm << " Out Edges : ";
                if (out_edges_.size() > 0)
                {
                    strm << out_edges_[0].Print();
                    for (size_t idx = 1; idx < out_edges_.size(); idx++)
                    {
                        strm << ", " << out_edges_[idx].Print();
                    }
                }
                return strm.str();
            }

            const NodeValueType& GetValueImmutable() const
            {
                return value_;
            }

            NodeValueType& GetValueMutable()
            {
                return value_;
            }

            void AddInEdge(const GraphEdge& new_in_edge)
            {
                in_edges_.push_back(new_in_edge);
            }

            void AddOutEdge(const GraphEdge& new_out_edge)
            {
                out_edges_.push_back(new_out_edge);
            }

            void AddEdgePair(const GraphEdge& new_in_edge, const GraphEdge& new_out_edge)
            {
                AddInEdge(new_in_edge);
                AddOutEdge(new_out_edge);
            }

            double GetDistance() const
            {
                return distance_;
            }

            void SetDistance(const double distance)
            {
                distance_ = distance;
            }

            const std::vector<GraphEdge>& GetInEdgesImmutable() const
            {
                return in_edges_;
            }

            std::vector<GraphEdge>& GetInEdgesMutable()
            {
                return in_edges_;
            }

            const std::vector<GraphEdge>& GetOutEdgesImmutable() const
            {
                return out_edges_;
            }

            std::vector<GraphEdge>& GetOutEdgesMutable()
            {
                return out_edges_;
            }

            void SetInEdges(const std::vector<GraphEdge>& new_in_edges)
            {
                in_edges_ = new_in_edges;
            }

            void SetOutEdges(const std::vector<GraphEdge>& new_out_edges)
            {
                out_edges_ = new_out_edges;
            }
    };

    template<typename NodeValueType, typename Allocator=std::allocator<NodeValueType>>
    class Graph
    {
        protected:

            std::vector<GraphNode<NodeValueType, Allocator>> nodes_;

        public:

            static uint64_t Serialize(const Graph<NodeValueType, Allocator>& graph, std::vector<uint8_t>& buffer, const std::function<uint64_t(const NodeValueType&, std::vector<uint8_t>&)>& value_serializer)
            {
                return graph.SerializeSelf(buffer, value_serializer);
            }

            static std::pair<Graph<NodeValueType, Allocator>, uint64_t> Deserialize(const std::vector<uint8_t>& buffer, const uint64_t current, const std::function<std::pair<NodeValueType, uint64_t>(const std::vector<uint8_t>&, const uint64_t)>& value_deserializer)
            {
                Graph<NodeValueType, Allocator> temp_graph;
                const uint64_t bytes_read = temp_graph.DeserializeSelf(buffer, current, value_deserializer);
                return std::make_pair(temp_graph, bytes_read);
            }

            Graph(const std::vector<GraphNode<NodeValueType, Allocator>>& nodes)
            {
                if (CheckGraphLinkage(nodes))
                {
                    nodes_ = nodes;
                }
                else
                {
                    throw std::invalid_argument("Invalid graph linkage");
                }
            }

            Graph(const size_t expected_size)
            {
                nodes_.reserve(expected_size);
            }

            Graph()
            {}

            uint64_t SerializeSelf(std::vector<uint8_t>& buffer, const std::function<uint64_t(const NodeValueType&, std::vector<uint8_t>&)>& value_serializer) const
            {
                const uint64_t start_buffer_size = buffer.size();
                std::function<uint64_t(const GraphNode<NodeValueType, Allocator>&, std::vector<uint8_t>&)> graph_state_serializer = std::bind(GraphNode<NodeValueType, Allocator>::Serialize, std::placeholders::_1, std::placeholders::_2, value_serializer);
                arc_helpers::SerializeVector<GraphNode<NodeValueType, Allocator>>(nodes_, buffer, graph_state_serializer);
                // Figure out how many bytes were written
                const uint64_t end_buffer_size = buffer.size();
                const uint64_t bytes_written = end_buffer_size - start_buffer_size;
                return bytes_written;
            }

            uint64_t DeserializeSelf(const std::vector<uint8_t>& buffer, const uint64_t current, const std::function<std::pair<NodeValueType, uint64_t>(const std::vector<uint8_t>&, const uint64_t)>& value_deserializer)
            {
                const std::function<std::pair<GraphNode<NodeValueType, Allocator>, uint64_t>(const std::vector<uint8_t>&, const uint64_t)> graph_state_deserializer = std::bind(GraphNode<NodeValueType, Allocator>::Deserialize, std::placeholders::_1, std::placeholders::_2, value_deserializer);
                const std::pair<std::vector<GraphNode<NodeValueType, Allocator>>, uint64_t> deserialized_nodes = arc_helpers::DeserializeVector<GraphNode<NodeValueType, Allocator>>(buffer, current, graph_state_deserializer);
                nodes_ = deserialized_nodes.first;
                return deserialized_nodes.second;
            }

            std::string Print() const
            {
                std::ostringstream strm;
                strm << "Graph - Nodes : ";
                if (nodes_.size() > 0)
                {
                    strm << nodes_[0].Print();
                    for (size_t idx = 1; idx < nodes_.size(); idx++)
                    {
                        strm << "\n" << nodes_[idx].Print();
                    }
                }
                return strm.str();
            }

            void ShrinkToFit()
            {
                nodes_.shrink_to_fit();
            }

            bool CheckGraphLinkage() const
            {
                return CheckGraphLinkage(GetNodesImmutable());
            }

            static bool CheckGraphLinkage(const Graph<NodeValueType, Allocator>& graph)
            {
                return CheckGraphLinkage(graph.GetNodesImmutable());
            }

            static bool CheckGraphLinkage(const std::vector<GraphNode<NodeValueType, Allocator>>& nodes)
            {
                // Go through every node and make sure the edges are valid
                for (size_t idx = 0; idx < nodes.size(); idx++)
                {
                    const GraphNode<NodeValueType, Allocator>& current_node = nodes[idx];
                    // Check the in edges first
                    const std::vector<GraphEdge>& in_edges = current_node.GetInEdgesImmutable();
                    for (size_t in_edge_idx = 0; in_edge_idx < in_edges.size(); in_edge_idx++)
                    {
                        const GraphEdge& current_edge = in_edges[in_edge_idx];
                        // Check from index to make sure it's in bounds
                        const int64_t from_index = current_edge.GetFromIndex();
                        if (from_index < 0 || from_index >= (int64_t)nodes.size())
                        {
                            return false;
                        }
                        // Check to index to make sure it matches our own index
                        const int64_t to_index = current_edge.GetToIndex();
                        if (to_index != (int64_t)idx)
                        {
                            return false;
                        }
                        // Check edge validity (edges to ourself are not allowed)
                        if (from_index == to_index)
                        {
                            return false;
                        }
                        // Check to make sure that the from index node is linked to us
                        const GraphNode<NodeValueType, Allocator>& from_node = nodes[(size_t)from_index];
                        const std::vector<GraphEdge>& from_node_out_edges = from_node.GetOutEdgesImmutable();
                        bool from_node_connection_valid = false;
                        // Make sure at least one out edge of the from index node corresponds to the current node
                        for (size_t from_node_out_edge_idx = 0; from_node_out_edge_idx < from_node_out_edges.size(); from_node_out_edge_idx++)
                        {
                            const GraphEdge& current_from_node_out_edge = from_node_out_edges[from_node_out_edge_idx];
                            if (current_from_node_out_edge.GetToIndex() == (int64_t)idx)
                            {
                                from_node_connection_valid = true;
                            }
                        }
                        if (from_node_connection_valid == false)
                        {
                            return false;
                        }
                    }
                    // Check the out edges second
                    const std::vector<GraphEdge>& out_edges = current_node.GetOutEdgesImmutable();
                    for (size_t out_edge_idx = 0; out_edge_idx < out_edges.size(); out_edge_idx++)
                    {
                        const GraphEdge& current_edge = out_edges[out_edge_idx];
                        // Check from index to make sure it matches our own index
                        const int64_t from_index = current_edge.GetFromIndex();
                        if (from_index != (int64_t)idx)
                        {
                            return false;
                        }
                        // Check to index to make sure it's in bounds
                        const int64_t to_index = current_edge.GetToIndex();
                        if (to_index < 0 || to_index >= (int64_t)nodes.size())
                        {
                            return false;
                        }
                        // Check edge validity (edges to ourself are not allowed)
                        if (from_index == to_index)
                        {
                            return false;
                        }
                        // Check to make sure that the to index node is linked to us
                        const GraphNode<NodeValueType, Allocator>& to_node = nodes[(size_t)to_index];
                        const std::vector<GraphEdge>& to_node_in_edges = to_node.GetInEdgesImmutable();
                        bool to_node_connection_valid = false;
                        // Make sure at least one in edge of the to index node corresponds to the current node
                        for (size_t to_node_in_edge_idx = 0; to_node_in_edge_idx < to_node_in_edges.size(); to_node_in_edge_idx++)
                        {
                            const GraphEdge& current_to_node_in_edge = to_node_in_edges[to_node_in_edge_idx];
                            if (current_to_node_in_edge.GetFromIndex() == (int64_t)idx)
                            {
                                to_node_connection_valid = true;
                            }
                        }
                        if (to_node_connection_valid == false)
                        {
                            return false;
                        }
                    }
                }
                return true;
            }

            const std::vector<GraphNode<NodeValueType, Allocator>>& GetNodesImmutable() const
            {
                return nodes_;
            }

            std::vector<GraphNode<NodeValueType, Allocator>>& GetNodesMutable()
            {
                return nodes_;
            }

            const GraphNode<NodeValueType, Allocator>& GetNodeImmutable(const int64_t index) const
            {
                assert(index >= 0);
                assert(index < (int64_t)nodes_.size());
                return nodes_[(size_t)index];
            }

            GraphNode<NodeValueType, Allocator>& GetNodeMutable(const int64_t index)
            {
                assert(index >= 0);
                assert(index < (int64_t)nodes_.size());
                return nodes_[(size_t)index];
            }

            int64_t AddNode(const GraphNode<NodeValueType, Allocator>& new_node)
            {
                nodes_.push_back(new_node);
                return (int64_t)(nodes_.size() - 1);
            }

            int64_t AddNode(const NodeValueType& new_value)
            {
                nodes_.push_back(GraphNode<NodeValueType, Allocator>(new_value));
                return (int64_t)(nodes_.size() - 1);
            }

            void AddEdgeBetweenNodes(const int64_t from_index, const int64_t to_index, const double edge_weight)
            {
                assert(from_index >= 0);
                assert(from_index < (int64_t)nodes_.size());
                assert(to_index >= 0);
                assert(to_index < (int64_t)nodes_.size());
                assert(from_index != to_index);
                const GraphEdge new_edge(from_index, to_index, edge_weight);
                GetNodeMutable(from_index).AddOutEdge(new_edge);
                GetNodeMutable(to_index).AddInEdge(new_edge);
            }

            void AddEdgesBetweenNodes(const int64_t first_index, const int64_t second_index, const double edge_weight)
            {
                assert(first_index >= 0);
                assert(first_index < (int64_t)nodes_.size());
                assert(second_index >= 0);
                assert(second_index < (int64_t)nodes_.size());
                assert(first_index != second_index);
                const GraphEdge first_edge(first_index, second_index, edge_weight);
                GetNodeMutable(first_index).AddOutEdge(first_edge);
                GetNodeMutable(second_index).AddInEdge(first_edge);
                const GraphEdge second_edge(second_index, first_index, edge_weight);
                GetNodeMutable(second_index).AddOutEdge(second_edge);
                GetNodeMutable(first_index).AddInEdge(second_edge);
            }
    };

    template<typename NodeValueType, typename Allocator=std::allocator<NodeValueType>>
    class SimpleDijkstrasAlgorithm
    {
        protected:

            class CompareIndexFn
            {
                public:

                    constexpr bool operator()(const std::pair<int64_t, double>& lhs, const std::pair<int64_t, double>& rhs) const
                    {
                        return lhs.second > rhs.second;
                    }
            };

            SimpleDijkstrasAlgorithm()
            {}

        public:
            typedef std::pair<Graph<NodeValueType, Allocator>, std::pair<std::vector<int64_t>, std::vector<double>>> DijkstrasResult;

            static DijkstrasResult PerformDijkstrasAlgorithm(const Graph<NodeValueType, Allocator>& graph, const int64_t start_index)
            {
                assert(start_index >= (int64_t)0);
                assert(start_index < (int64_t)graph.GetNodesImmutable().size());
                Graph<NodeValueType, Allocator> working_copy = graph;
                // Setup
                std::vector<int64_t> previous_index_map(working_copy.GetNodesImmutable().size(), -1);
                std::vector<double> distances(working_copy.GetNodesImmutable().size(), std::numeric_limits<double>::infinity());
                std::priority_queue<std::pair<int64_t, double>, std::vector<std::pair<int64_t, double>>, CompareIndexFn> queue;
                std::unordered_map<int64_t, uint32_t> explored(graph.GetNodesImmutable().size());
                for (size_t idx = 0; idx < working_copy.GetNodesImmutable().size(); idx++)
                {
                    working_copy.GetNodeMutable((int64_t)idx).SetDistance(std::numeric_limits<double>::infinity());
                    queue.push(std::make_pair((int64_t)idx, std::numeric_limits<double>::infinity()));
                }
                working_copy.GetNodeMutable(start_index).SetDistance(0.0);
                previous_index_map[(size_t)start_index] = start_index;
                distances[(size_t)start_index] = 0.0;
                queue.push(std::make_pair(start_index, 0.0));
                while (queue.size() > 0)
                {
                    const std::pair<int64_t, double> top_node = queue.top();
                    const int64_t& top_node_index = top_node.first;
                    const double& top_node_distance = top_node.second;
                    queue.pop();
                    if (explored[top_node.first] > 0)
                    {
                        // We've already been here
                        continue;
                    }
                    else
                    {
                        // Note that we've been here
                        explored[top_node.first] = 1;
                        // Get our neighbors
                        const std::vector<GraphEdge>& neighbor_edges = working_copy.GetNodeImmutable(top_node_index).GetInEdgesImmutable();
                        // Go through our neighbors
                        for (size_t neighbor_idx = 0; neighbor_idx < neighbor_edges.size(); neighbor_idx++)
                        {
                            const int64_t neighbor_index = neighbor_edges[neighbor_idx].GetFromIndex();
                            const double neighbor_edge_weight = neighbor_edges[neighbor_idx].GetWeight();
                            const double new_neighbor_distance = top_node_distance + neighbor_edge_weight;
                            // Check against the neighbor
                            const double stored_neighbor_distance = working_copy.GetNodeImmutable(neighbor_index).GetDistance();
                            if (new_neighbor_distance < stored_neighbor_distance)
                            {
                                // We've found a better way to get to this node
                                // Check if it's already been explored
                                if (explored[neighbor_index] > 0)
                                {
                                    // If it's already been explored, we just update it in place
                                    working_copy.GetNodeMutable(neighbor_index).SetDistance(new_neighbor_distance);
                                }
                                else
                                {
                                    // If it hasn't been explored, we need to update it and add it to the queue
                                    working_copy.GetNodeMutable(neighbor_index).SetDistance(new_neighbor_distance);
                                    queue.push(std::make_pair(neighbor_index, new_neighbor_distance));
                                }
                                // Update that we're the best previous node
                                previous_index_map[(size_t)neighbor_index] = top_node_index;
                                distances[(size_t)neighbor_index] = new_neighbor_distance;
                            }
                            else
                            {
                                // Do nothing
                                continue;
                            }
                        }
                    }
                }
                return std::make_pair(working_copy, std::make_pair(previous_index_map, distances));
            }

            // These functions have not been tested.  Use with care.
            static uint64_t SerializeDijstrasResult(const DijkstrasResult& result, std::vector<uint8_t>& buffer, const std::function<uint64_t(const NodeValueType&, std::vector<uint8_t>&)>& value_serializer)
            {
                const uint64_t start_buffer_size = buffer.size();
                // Serialize the graph
                result.first.SerializeSelf(buffer, value_serializer);
                // Serialize the previous index
                SerializeVector(result.second.first, std::bind(arc_helpers::SerializeFixedSizePOD<uint64_t>, std::placeholders::_1, std::placeholders::_2));
                // Serialze the distances
                SerializeVector(result.second.second, std::bind(arc_helpers::SerializeFixedSizePOD<uint64_t>, std::placeholders::_1, std::placeholders::_2));
                // Figure out how many bytes were written
                const uint64_t end_buffer_size = buffer.size();
                const uint64_t bytes_written = end_buffer_size - start_buffer_size;
                return bytes_written;
            }

            // These functions have not been tested.  Use with care.
            static std::pair<DijkstrasResult, uint64_t> DijstrasResult(const std::vector<uint8_t>& buffer, const uint64_t current, const std::function<std::pair<NodeValueType, uint64_t>(const std::vector<uint8_t>&, const uint64_t)>& value_deserializer)
            {
                assert(current < buffer.size());
                uint64_t current_position = current;
                // Deserialize the graph itself
                std::pair<DijkstrasResult, uint64_t> deserialized;
                const std::pair<Graph<NodeValueType, Allocator>, uint64_t> graph_deserialized = Graph<NodeValueType, Allocator>::Deserialize(buffer, current_position, value_deserializer);
                deserialized.first.first = graph_deserialized.first;
                current_position += graph_deserialized.second;
                // Deserialize the previous index
                const std::pair<std::vector<int64_t>, uint64_t> prev_index_deserialized = arc_helpers::DeserializeVector<int64_t>(buffer, current_position, std::bind(arc_helpers::DeserializeFixedSizePOD<uint64_t>, std::placeholders::_1, std::placeholders::_2));
                deserialized.first.second.first = prev_index_deserialized.first;
                current_position += prev_index_deserialized.second;
                // Deserialize the distances
                const std::pair<std::vector<double>, uint64_t> distance_deserialized = arc_helpers::DeserializeVector<double>(buffer, current_position, std::bind(arc_helpers::DeserializeFixedSizePOD<double>, std::placeholders::_1, std::placeholders::_2));
                deserialized.first.second.second = distance_deserialized.first;
                current_position += distance_deserialized.second;
                // Figure out how many bytes were read
                deserialized.second = current_position - current;
                return deserialized;
            }
    };

}

#endif // DIJKSTRAS_HPP
