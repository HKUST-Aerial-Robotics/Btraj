#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <random>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <visualization_msgs/Marker.h>
#include <arc_utilities/voxel_grid.hpp>
#include <arc_utilities/arc_helpers.hpp>
#include <arc_utilities/eigen_helpers.hpp>
#include <arc_utilities/pretty_print.hpp>
#include <arc_utilities/simple_kmeans_clustering.hpp>
#include <sdf_tools/sdf.hpp>
#include <sdf_tools/TaggedObjectCollisionMap.h>

#ifndef TAGGED_OBJECT_COLLISION_MAP_HPP
#define TAGGED_OBJECT_COLLISION_MAP_HPP

#define ENABLE_UNORDERED_MAP_SIZE_HINTS

namespace sdf_tools
{
    struct TAGGED_OBJECT_COLLISION_CELL
    {
        float occupancy;
        uint32_t component;
        uint32_t object_id;
        uint32_t convex_segment;

        TAGGED_OBJECT_COLLISION_CELL() : occupancy(0.0), component(0u), object_id(0u), convex_segment(0u) {}

        TAGGED_OBJECT_COLLISION_CELL(const float in_occupancy, const uint32_t in_object_id) : occupancy(in_occupancy), component(0u), object_id(in_object_id), convex_segment(0u) {}

        TAGGED_OBJECT_COLLISION_CELL(const float in_occupancy, const uint32_t in_object_id, const uint32_t in_component, const uint32_t in_convex_segment) : occupancy(in_occupancy), component(in_component), object_id(in_object_id), convex_segment(in_convex_segment) {}

        bool SharesConvexSegment(const TAGGED_OBJECT_COLLISION_CELL& other) const
        {
            if ((convex_segment & other.convex_segment) > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        std::vector<uint32_t> GetListOfConvexSegments() const
        {
            uint32_t temp_convex_segment = convex_segment;
            std::vector<uint32_t> convex_segments;
            for (uint32_t segment = 1; segment <= 32; segment++)
            {
                if ((temp_convex_segment & 0x00000001) == 1)
                {
                    convex_segments.push_back(segment);
                }
                temp_convex_segment = temp_convex_segment >> 1;
            }
            return convex_segments;
        }

        bool IsPartOfConvexSegment(const uint32_t segment) const
        {
            assert(segment >= 1);
            assert(segment <= 32);
            const uint32_t mask = arc_helpers::SetBit(0u, segment - 1, true);
            if ((mask & convex_segment) > 0)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        void AddToConvexSegment(const uint32_t segment)
        {
            assert(segment >= 1);
            assert(segment <= 32);
            convex_segment = arc_helpers::SetBit(convex_segment, segment - 1, true);
        }

        void RemoveFromConvexSegment(const uint32_t segment)
        {
            assert(segment >= 1);
            assert(segment <= 32);
            convex_segment = arc_helpers::SetBit(convex_segment, segment - 1, false);
        }
    };

    inline std::vector<uint8_t> TaggedObjectCollisionCellToBinary(const TAGGED_OBJECT_COLLISION_CELL& value)
    {
        std::vector<uint8_t> binary(sizeof(TAGGED_OBJECT_COLLISION_CELL));
        memcpy(&binary.front(), &value, sizeof(TAGGED_OBJECT_COLLISION_CELL));
        return binary;
    }

    inline TAGGED_OBJECT_COLLISION_CELL TaggedObjectCollisionCellFromBinary(const std::vector<uint8_t>& binary)
    {
        if (binary.size() != sizeof(TAGGED_OBJECT_COLLISION_CELL))
        {
            std::cerr << "Binary value is not " << sizeof(TAGGED_OBJECT_COLLISION_CELL) << " bytes" << std::endl;
            return TAGGED_OBJECT_COLLISION_CELL();
        }
        else
        {
            TAGGED_OBJECT_COLLISION_CELL loaded;
            memcpy(&loaded, &binary.front(), sizeof(TAGGED_OBJECT_COLLISION_CELL));
            return loaded;
        }
    }

    class TaggedObjectCollisionMapGrid
    {
    protected:

        inline static std_msgs::ColorRGBA GenerateComponentColor(const uint32_t component, const float alpha=1.0f)
        {
            return arc_helpers::GenerateUniqueColor<std_msgs::ColorRGBA>(component, alpha);
        }

        inline bool IsSurfaceIndex(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            // First, we make sure that indices are within bounds
            // Out of bounds indices are NOT surface cells
            if (x_index < 0 || y_index < 0 || z_index < 0 || x_index >= GetNumXCells() || y_index >= GetNumYCells() || z_index >= GetNumZCells())
            {
                return false;
            }
            // Edge indices are automatically surface cells
            if (x_index == 0 || y_index == 0 || z_index == 0 || x_index == (GetNumXCells() - 1) || y_index == (GetNumYCells() - 1) || z_index == (GetNumZCells()))
            {
                return true;
            }
            // If the cell is inside the grid, we check the neighbors
            // Note that we must check all 26 neighbors
            uint32_t our_component = GetImmutable(x_index, y_index, z_index).first.component;
            // Check neighbor 1
            if (our_component != GetImmutable(x_index, y_index, z_index - 1).first.component)
            {
                return true;
            }
            // Check neighbor 2
            else if (our_component != GetImmutable(x_index, y_index, z_index + 1).first.component)
            {
                return true;
            }
            // Check neighbor 3
            else if (our_component != GetImmutable(x_index, y_index - 1, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 4
            else if (our_component != GetImmutable(x_index, y_index + 1, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 5
            else if (our_component != GetImmutable(x_index - 1, y_index, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 6
            else if (our_component != GetImmutable(x_index + 1, y_index, z_index).first.component)
            {
                return true;
            }
            // If none of the faces are exposed, it's not a surface voxel
            return false;
        }

        struct bucket_cell
        {
            uint32_t location[3];
            uint32_t closest_point[3];
            double distance_square;
            int32_t update_direction;
        };

        typedef VoxelGrid::VoxelGrid<bucket_cell> DistanceField;

        inline DistanceField BuildDistanceField(const std::vector<VoxelGrid::GRID_INDEX>& points) const
        {
            // Make the DistanceField container
            bucket_cell default_cell;
            default_cell.distance_square = INFINITY;
            DistanceField distance_field(GetOriginTransform(), GetResolution(), GetXSize(), GetYSize(), GetZSize(), default_cell);
            // Compute maximum distance square
            long max_distance_square = (distance_field.GetNumXCells() * distance_field.GetNumXCells()) + (distance_field.GetNumYCells() * distance_field.GetNumYCells()) + (distance_field.GetNumZCells() * distance_field.GetNumZCells());
            // Make bucket queue
            std::vector<std::vector<bucket_cell>> bucket_queue(max_distance_square + 1);
            bucket_queue[0].reserve(points.size());
            // Set initial update direction
            int initial_update_direction = GetDirectionNumber(0, 0, 0);
            // Mark all points with distance zero and add to the bucket queue
            for (size_t index = 0; index < points.size(); index++)
            {
                const VoxelGrid::GRID_INDEX& current_index = points[index];
                std::pair<bucket_cell&, bool> query = distance_field.GetMutable(current_index);
                if (query.second)
                {
                    query.first.location[0] = current_index.x;
                    query.first.location[1] = current_index.y;
                    query.first.location[2] = current_index.z;
                    query.first.closest_point[0] = current_index.x;
                    query.first.closest_point[1] = current_index.y;
                    query.first.closest_point[2] = current_index.z;
                    query.first.distance_square = 0.0;
                    query.first.update_direction = initial_update_direction;
                    bucket_queue[0].push_back(query.first);
                }
                // If the point is outside the bounds of the SDF, skip
                else
                {
                    continue;
                }
            }
            // Process the bucket queue
            std::vector<std::vector<std::vector<std::vector<int>>>> neighborhoods = MakeNeighborhoods();
            for (size_t bq_idx = 0; bq_idx < bucket_queue.size(); bq_idx++)
            {
                std::vector<bucket_cell>::iterator queue_itr = bucket_queue[bq_idx].begin();
                while (queue_itr != bucket_queue[bq_idx].end())
                {
                    // Get the current location
                    bucket_cell& cur_cell = *queue_itr;
                    double x = cur_cell.location[0];
                    double y = cur_cell.location[1];
                    double z = cur_cell.location[2];
                    // Pick the update direction
                    int D = bq_idx;
                    if (D > 1)
                    {
                        D = 1;
                    }
                    // Make sure the update direction is valid
                    if (cur_cell.update_direction < 0 || cur_cell.update_direction > 26)
                    {
                        ++queue_itr;
                        continue;
                    }
                    // Get the current neighborhood list
                    std::vector<std::vector<int>>& neighborhood = neighborhoods[D][cur_cell.update_direction];
                    // Update the distance from the neighboring cells
                    for (size_t nh_idx = 0; nh_idx < neighborhood.size(); nh_idx++)
                    {
                        // Get the direction to check
                        int dx = neighborhood[nh_idx][0];
                        int dy = neighborhood[nh_idx][1];
                        int dz = neighborhood[nh_idx][2];
                        int nx = x + dx;
                        int ny = y + dy;
                        int nz = z + dz;
                        std::pair<bucket_cell&, bool> neighbor_query = distance_field.GetMutable((int64_t)nx, (int64_t)ny, (int64_t)nz);
                        if (!neighbor_query.second)
                        {
                            // "Neighbor" is outside the bounds of the SDF
                            continue;
                        }
                        // Update the neighbor's distance based on the current
                        int new_distance_square = ComputeDistanceSquared(nx, ny, nz, cur_cell.closest_point[0], cur_cell.closest_point[1], cur_cell.closest_point[2]);
                        if (new_distance_square > max_distance_square)
                        {
                            // Skip these cases
                            continue;
                        }
                        if (new_distance_square < neighbor_query.first.distance_square)
                        {
                            // If the distance is better, time to update the neighbor
                            neighbor_query.first.distance_square = new_distance_square;
                            neighbor_query.first.closest_point[0] = cur_cell.closest_point[0];
                            neighbor_query.first.closest_point[1] = cur_cell.closest_point[1];
                            neighbor_query.first.closest_point[2] = cur_cell.closest_point[2];
                            neighbor_query.first.location[0] = nx;
                            neighbor_query.first.location[1] = ny;
                            neighbor_query.first.location[2] = nz;
                            neighbor_query.first.update_direction = GetDirectionNumber(dx, dy, dz);
                            // Add the neighbor into the bucket queue
                            bucket_queue[new_distance_square].push_back(neighbor_query.first);
                        }
                    }
                    // Increment the queue iterator
                    ++queue_itr;
                }
                // Clear the current queue now that we're done with it
                bucket_queue[bq_idx].clear();
            }
            return distance_field;
        }

        inline std::vector<std::vector<std::vector<std::vector<int>>>> MakeNeighborhoods()  const
        {
            std::vector<std::vector<std::vector<std::vector<int>>>> neighborhoods;
            neighborhoods.resize(2);
            for (size_t n = 0; n < neighborhoods.size(); n++)
            {
                neighborhoods[n].resize(27);
                // Loop through the source directions
                for (int dx = -1; dx <= 1; dx++)
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        for (int dz = -1; dz <= 1; dz++)
                        {
                            int direction_number = GetDirectionNumber(dx, dy, dz);
                            // Loop through the target directions
                            for (int tdx = -1; tdx <= 1; tdx++)
                            {
                                for (int tdy = -1; tdy <= 1; tdy++)
                                {
                                    for (int tdz = -1; tdz <= 1; tdz++)
                                    {
                                        if (tdx == 0 && tdy == 0 && tdz == 0)
                                        {
                                            continue;
                                        }
                                        if (n >= 1)
                                        {
                                            if ((abs(tdx) + abs(tdy) + abs(tdz)) != 1)
                                            {
                                                continue;
                                            }
                                            if ((dx * tdx) < 0 || (dy * tdy) < 0 || (dz * tdz) < 0)
                                            {
                                                continue;
                                            }
                                        }
                                        std::vector<int> new_point;
                                        new_point.resize(3);
                                        new_point[0] = tdx;
                                        new_point[1] = tdy;
                                        new_point[2] = tdz;
                                        neighborhoods[n][direction_number].push_back(new_point);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            return neighborhoods;
        }

        inline int GetDirectionNumber(const int dx, const int dy, const int dz) const
        {
            return ((dx + 1) * 9) + ((dy + 1) * 3) + (dz + 1);
        }

        inline double ComputeDistanceSquared(const int32_t x1, const int32_t y1, const int32_t z1, const int32_t x2, const int32_t y2, const int32_t z2) const
        {
            int32_t dx = x1 - x2;
            int32_t dy = y1 - y2;
            int32_t dz = z1 - z2;
            return double((dx * dx) + (dy * dy) + (dz * dz));
        }

        VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> collision_field_;
        uint32_t number_of_components_;
        std::string frame_;
        bool initialized_;
        bool components_valid_;
        bool convex_segments_valid_;

        std::vector<uint8_t> PackBinaryRepresentation(const std::vector<TAGGED_OBJECT_COLLISION_CELL>& raw) const;

        std::vector<TAGGED_OBJECT_COLLISION_CELL> UnpackBinaryRepresentation(const std::vector<uint8_t>& packed) const;

        int64_t MarkConnectedComponent(const int64_t x_index, const int64_t y_index, const int64_t z_index, const uint32_t connected_component)
        {
            // Make the working queue
            std::list<VoxelGrid::GRID_INDEX> working_queue;
            // Make a hash table to store queued indices (so we don't repeat work)
            // Let's provide an hint at the size of hashmap we'll need, since this will reduce the need to resize & rehash
            // We're going to assume that connected components, in general, will take ~1/16 of the grid in size
            // which means, with 2 cells/hash bucket, we'll initialize to grid size/32
        #ifdef ENABLE_UNORDERED_MAP_SIZE_HINTS
            size_t queued_hashtable_size_hint = collision_field_.GetRawData().size() / 32;
            std::unordered_map<VoxelGrid::GRID_INDEX, int8_t> queued_hashtable(queued_hashtable_size_hint);
        #else
            std::unordered_map<VoxelGrid::GRID_INDEX, int8_t> queued_hashtable;
        #endif
            // Add the starting index
            VoxelGrid::GRID_INDEX start_index(x_index, y_index, z_index);
            // Enqueue it
            working_queue.push_back(start_index);
            queued_hashtable[start_index] = 1;
            // Work
            int64_t marked_cells = 0;
            while (working_queue.size() > 0)
            {
                // Get an item off the queue to work with
                VoxelGrid::GRID_INDEX current_index = working_queue.front();
                working_queue.pop_front();
                // Get the current value
                TAGGED_OBJECT_COLLISION_CELL& current_value = GetMutable(current_index.x, current_index.y, current_index.z).first;
                // Mark the connected component
                current_value.component = connected_component;
        //        // Update the grid
        //        Set(current_index.x, current_index.y, current_index.z, current_value);
                // Go through the possible neighbors and enqueue as needed
                // Since there are only six cases (voxels must share a face to be considered connected), we handle each explicitly
                // Case 1
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1_neighbor = GetImmutable(current_index.x - 1, current_index.y, current_index.z);
                if (xm1_neighbor.second && (current_value.occupancy == xm1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x - 1, current_index.y, current_index.z);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
                // Case 2
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> ym1_neighbor = GetImmutable(current_index.x, current_index.y - 1, current_index.z);
                if (ym1_neighbor.second && (current_value.occupancy == ym1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x, current_index.y - 1, current_index.z);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
                // Case 3
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> zm1_neighbor = GetImmutable(current_index.x, current_index.y, current_index.z - 1);
                if (zm1_neighbor.second && (current_value.occupancy == zm1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x, current_index.y, current_index.z - 1);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
                // Case 4
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1_neighbor = GetImmutable(current_index.x + 1, current_index.y, current_index.z);
                if (xp1_neighbor.second && (current_value.occupancy == xp1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x + 1, current_index.y, current_index.z);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
                // Case 5
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> yp1_neighbor = GetImmutable(current_index.x, current_index.y + 1, current_index.z);
                if (yp1_neighbor.second && (current_value.occupancy == yp1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x, current_index.y + 1, current_index.z);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
                // Case 6
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> zp1_neighbor = GetImmutable(current_index.x, current_index.y, current_index.z + 1);
                if (zp1_neighbor.second && (current_value.occupancy == zp1_neighbor.first.occupancy))
                {
                    VoxelGrid::GRID_INDEX neighbor_index(current_index.x, current_index.y, current_index.z + 1);
                    if (queued_hashtable[neighbor_index] <= 0)
                    {
                        queued_hashtable[neighbor_index] = 1;
                        working_queue.push_back(neighbor_index);
                    }
                }
            }
            return marked_cells;
        }

        std::vector<VoxelGrid::GRID_INDEX> CheckIfConvex(const VoxelGrid::GRID_INDEX& candidate_index, std::unordered_map<VoxelGrid::GRID_INDEX, int8_t>& explored_indices, const VoxelGrid::VoxelGrid<std::vector<uint32_t>>& region_grid, const uint32_t current_convex_region) const;

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        inline TaggedObjectCollisionMapGrid(const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const TAGGED_OBJECT_COLLISION_CELL& default_value, const TAGGED_OBJECT_COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> new_field(resolution, x_size, y_size, z_size, default_value, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
            convex_segments_valid_ = false;
        }

        inline TaggedObjectCollisionMapGrid(const Eigen::Affine3d& origin_transform, const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const TAGGED_OBJECT_COLLISION_CELL& default_value, const TAGGED_OBJECT_COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> new_field(origin_transform, resolution, x_size, y_size, z_size, default_value, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
            convex_segments_valid_ = false;
        }

        inline TaggedObjectCollisionMapGrid(const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const TAGGED_OBJECT_COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> new_field(resolution, x_size, y_size, z_size, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
            convex_segments_valid_ = false;
        }

        inline TaggedObjectCollisionMapGrid(const Eigen::Affine3d& origin_transform, const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const TAGGED_OBJECT_COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> new_field(origin_transform, resolution, x_size, y_size, z_size, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
            convex_segments_valid_ = false;
        }

        inline TaggedObjectCollisionMapGrid() : number_of_components_(0), initialized_(false), components_valid_(false), convex_segments_valid_(false) {}

        inline bool IsInitialized() const
        {
            return initialized_;
        }

        inline bool AreComponentsValid() const
        {
            return components_valid_;
        }

        inline bool AreConvexSegmentsValid() const
        {
            return convex_segments_valid_;
        }

        inline std::pair<bool, bool> CheckIfCandidateCorner3d(const Eigen::Vector3d& location) const
        {
            const std::vector<int64_t> indices = collision_field_.LocationToGridIndex3d(location);
            if (indices.size() == 3)
            {
                return CheckIfCandidateCorner(indices[0], indices[1], indices[2]);
            }
            else
            {
                return std::pair<bool, bool>(false, false);
            }
        }

        inline std::pair<bool, bool> CheckIfCandidateCorner4d(const Eigen::Vector4d& location) const
        {
            const std::vector<int64_t> indices = collision_field_.LocationToGridIndex4d(location);
            if (indices.size() == 3)
            {
                return CheckIfCandidateCorner(indices[0], indices[1], indices[2]);
            }
            else
            {
                return std::pair<bool, bool>(false, false);
            }
        }

        inline std::pair<bool, bool> CheckIfCandidateCorner(const double x, const double y, const double z) const
        {
            const Eigen::Vector4d location(x, y, z, 1.0);
            return CheckIfCandidateCorner4d(location);
        }

        inline std::pair<bool, bool> CheckIfCandidateCorner(const VoxelGrid::GRID_INDEX& index) const
        {
            return CheckIfCandidateCorner(index.x, index.y, index.z);
        }

        inline std::pair<bool, bool> CheckIfCandidateCorner(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            assert(components_valid_);
            assert(convex_segments_valid_);
            const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> current_cell = GetImmutable(x_index, y_index, z_index);
            if (current_cell.second)
            {
                // Grab the six neighbors ONLY if they belong to a different component
                uint32_t different_neighbors = 0u;
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1yz_cell = GetImmutable(x_index - 1, y_index, z_index);
                if (xm1yz_cell.second && (xm1yz_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1yz_cell = GetImmutable(x_index + 1, y_index, z_index);
                if (xp1yz_cell.second && (xp1yz_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xym1z_cell = GetImmutable(x_index, y_index - 1, z_index);
                if (xym1z_cell.second && (xym1z_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyp1z_cell = GetImmutable(x_index, y_index + 1, z_index);
                if (xyp1z_cell.second && (xyp1z_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyzm1_cell = GetImmutable(x_index, y_index, z_index - 1);
                if (xyzm1_cell.second && (xyzm1_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                const std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyzp1_cell = GetImmutable(x_index, y_index, z_index + 1);
                if (xyzp1_cell.second && (xyzp1_cell.first.component != current_cell.first.component))
                {
                    different_neighbors++;
                }
                // We now have between zero and six neighbors to work with
                if (different_neighbors <= 1u)
                {
                    // If there is one or fewer neighbors to work with, we are clearly not a corner
                    return std::pair<bool, bool>(false, true);
                }
                else
                {
                    // If there are 2 or more neighbors to work with, we are a candidate corner
                    return std::pair<bool, bool>(true, true);
                }
            }
            else
            {
                // Not in the grid
                return std::pair<bool, bool>(false, false);
            }
        }

        inline std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> GetImmutable3d(const Eigen::Vector3d& location) const
        {
            return collision_field_.GetImmutable3d(location);
        }

        inline std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> GetImmutable4d(const Eigen::Vector4d& location) const
        {
            return collision_field_.GetImmutable4d(location);
        }

        inline std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> GetImmutable(const double x, const double y, const double z) const
        {
            return collision_field_.GetImmutable(x, y, z);
        }

        inline std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> GetImmutable(const VoxelGrid::GRID_INDEX& index) const
        {
            return collision_field_.GetImmutable(index);
        }

        inline std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> GetImmutable(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return collision_field_.GetImmutable(x_index, y_index, z_index);
        }

        inline std::pair<TAGGED_OBJECT_COLLISION_CELL&, bool> GetMutable3d(const Eigen::Vector3d& location)
        {
            return collision_field_.GetMutable3d(location);
        }

        inline std::pair<TAGGED_OBJECT_COLLISION_CELL&, bool> GetMutable4d(const Eigen::Vector4d& location)
        {
            return collision_field_.GetMutable4d(location);
        }

        inline std::pair<TAGGED_OBJECT_COLLISION_CELL&, bool> GetMutable(const double x, const double y, const double z)
        {
            return collision_field_.GetMutable(x, y, z);
        }

        inline std::pair<TAGGED_OBJECT_COLLISION_CELL&, bool> GetMutable(const VoxelGrid::GRID_INDEX& index)
        {
            return collision_field_.GetMutable(index);
        }

        inline std::pair<TAGGED_OBJECT_COLLISION_CELL&, bool> GetMutable(const int64_t x_index, const int64_t y_index, const int64_t z_index)
        {
            return collision_field_.GetMutable(x_index, y_index, z_index);
        }

        inline bool Set(const double x, const double y, const double z, const TAGGED_OBJECT_COLLISION_CELL& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(x, y, z, value);
        }

        inline bool Set3d(const Eigen::Vector3d& location, const TAGGED_OBJECT_COLLISION_CELL& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue3d(location, value);
        }

        inline bool Set4d(const Eigen::Vector4d& location, const TAGGED_OBJECT_COLLISION_CELL& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue4d(location, value);
        }

        inline bool Set(const int64_t x_index, const int64_t y_index, const int64_t z_index, const TAGGED_OBJECT_COLLISION_CELL& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(x_index, y_index, z_index, value);
        }

        inline bool Set(const VoxelGrid::GRID_INDEX& index, const TAGGED_OBJECT_COLLISION_CELL& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(index, value);
        }

        inline bool Set(const double x, const double y, const double z, TAGGED_OBJECT_COLLISION_CELL&& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(x, y, z, value);
        }

        inline bool Set3d(const Eigen::Vector3d& location, TAGGED_OBJECT_COLLISION_CELL&& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue3d(location, value);
        }

        inline bool Set4d(const Eigen::Vector4d& location, TAGGED_OBJECT_COLLISION_CELL&& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue4d(location, value);
        }

        inline bool Set(const int64_t x_index, const int64_t y_index, const int64_t z_index, TAGGED_OBJECT_COLLISION_CELL&& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(x_index, y_index, z_index, value);
        }

        inline bool Set(const VoxelGrid::GRID_INDEX& index, TAGGED_OBJECT_COLLISION_CELL&& value)
        {
            components_valid_ = false;
            convex_segments_valid_ = false;
            return collision_field_.SetValue(index, value);
        }

        inline double GetXSize() const
        {
            return collision_field_.GetXSize();
        }

        inline double GetYSize() const
        {
            return collision_field_.GetYSize();
        }

        inline double GetZSize() const
        {
            return collision_field_.GetZSize();
        }

        inline double GetResolution() const
        {
            return collision_field_.GetCellSizes()[0];
        }

        inline TAGGED_OBJECT_COLLISION_CELL GetDefaultValue() const
        {
            return collision_field_.GetDefaultValue();
        }

        inline TAGGED_OBJECT_COLLISION_CELL GetOOBValue() const
        {
            return collision_field_.GetOOBValue();
        }

        inline int64_t GetNumXCells() const
        {
            return collision_field_.GetNumXCells();
        }

        inline int64_t GetNumYCells() const
        {
            return collision_field_.GetNumYCells();
        }

        inline int64_t GetNumZCells() const
        {
            return collision_field_.GetNumZCells();
        }

        inline const Eigen::Affine3d& GetOriginTransform() const
        {
            return collision_field_.GetOriginTransform();
        }

        inline const Eigen::Affine3d& GetInverseOriginTransform() const
        {
            return collision_field_.GetInverseOriginTransform();
        }

        inline std::string GetFrame() const
        {
            return frame_;
        }

        inline std::pair<uint32_t, bool> GetNumConnectedComponents() const
        {
            return std::pair<uint32_t, bool>(number_of_components_, components_valid_);
        }

        inline std::vector<int64_t> LocationToGridIndex3d(const Eigen::Vector3d& location) const
        {
            return collision_field_.LocationToGridIndex3d(location);
        }

        inline std::vector<int64_t> LocationToGridIndex4d(const Eigen::Vector4d& location) const
        {
            return collision_field_.LocationToGridIndex4d(location);
        }

        inline std::vector<int64_t> LocationToGridIndex(const double x, const double y, const double z) const
        {
            return collision_field_.LocationToGridIndex(x, y, z);
        }

        inline std::vector<double> GridIndexToLocation(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return collision_field_.GridIndexToLocation(x_index, y_index, z_index);
        }

        bool SaveToFile(const std::string& filepath) const;

        bool LoadFromFile(const std::string &filepath);

        sdf_tools::TaggedObjectCollisionMap GetMessageRepresentation() const;

        bool LoadFromMessageRepresentation(const sdf_tools::TaggedObjectCollisionMap& message);

        TaggedObjectCollisionMapGrid Resample(const double new_resolution) const
        {
            TaggedObjectCollisionMapGrid resampled(GetOriginTransform(), GetFrame(), new_resolution, GetXSize(), GetYSize(), GetZSize(), GetOOBValue());
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                        const std::vector<double> raw_current_cell_location = GridIndexToLocation(x_index, y_index, z_index);
                        assert(raw_current_cell_location.size() == 3);
                        const Eigen::Vector4d current_cell_location(raw_current_cell_location[0], raw_current_cell_location[1], raw_current_cell_location[2], 1.0);
                        resampled.Set4d(current_cell_location, current_cell);
                    }
                }
            }
            return resampled;
        }

        uint32_t UpdateConnectedComponents()
        {
            // If the connected components are already valid, skip computing them again
            if (components_valid_)
            {
                return number_of_components_;
            }
            components_valid_ = false;
            // Reset components first
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        GetMutable(x_index, y_index, z_index).first.component = 0;
                    }
                }
            }
            // Mark the components
            int64_t total_cells = GetNumXCells() * GetNumYCells() * GetNumZCells();
            int64_t marked_cells = 0;
            uint32_t connected_components = 0;
            // Sweep through the grid
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        // Check if the cell has already been marked, if so, ignore
                        if (GetImmutable(x_index, y_index, z_index).first.component > 0)
                        {
                            continue;
                        }
                        // Start marking a new connected component
                        connected_components++;
                        int64_t cells_marked = MarkConnectedComponent(x_index, y_index, z_index, connected_components);
                        marked_cells += cells_marked;
                        // Short-circuit if we've marked everything
                        if (marked_cells == total_cells)
                        {
                            number_of_components_ = connected_components;
                            components_valid_ = true;
                            return connected_components;
                        }
                    }
                }
            }
            number_of_components_ = connected_components;
            components_valid_ = true;
            return connected_components;
        }

        std::map<uint32_t, std::pair<int32_t, int32_t>> ComputeComponentTopology(const bool ignore_empty_components, const bool recompute_connected_components, const bool verbose)
        {
            // Recompute the connected components if need be
            if (recompute_connected_components)
            {
                UpdateConnectedComponents();
            }
            // Extract the surfaces of each connected component
            std::map<uint32_t, std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>> component_surfaces = ExtractComponentSurfaces(ignore_empty_components);
            // Compute the number of holes in each surface
            std::map<uint32_t, std::pair<int32_t, int32_t>> component_holes;
            std::map<uint32_t, std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>>::iterator component_surfaces_itr;
            for (component_surfaces_itr = component_surfaces.begin(); component_surfaces_itr != component_surfaces.end(); ++component_surfaces_itr)
            {
                uint32_t component_number = component_surfaces_itr->first;
                std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& component_surface = component_surfaces_itr->second;
                std::pair<int32_t, int32_t> number_of_holes_and_voids = ComputeHolesInSurface(component_number, component_surface, verbose);
                component_holes[component_number] = number_of_holes_and_voids;
            }
            return component_holes;
        }

        std::map<uint32_t, std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>> ExtractComponentSurfaces(const bool ignore_empty_components) const
        {
            std::map<uint32_t, std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>> component_surfaces;
            // Loop through the grid and extract surface cells for each component
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                        if (ignore_empty_components)
                        {
                            if (current_cell.occupancy > 0.5)
                            {
                                VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
                                if (IsSurfaceIndex(x_index, y_index, z_index))
                                {
                                    component_surfaces[current_cell.component][current_index] = 1;
                                }
                            }
                        }
                        else
                        {
                            VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
                            if (IsSurfaceIndex(x_index, y_index, z_index))
                            {
                                component_surfaces[current_cell.component][current_index] = 1;
                            }
                        }
                    }
                }
            }
            return component_surfaces;
        }

        /* Extracts the active indices from a surface map as a vector, which is useful in contexts where a 1-dimensional index into the surface is needed
         */
        std::vector<VoxelGrid::GRID_INDEX> ExtractStaticSurface(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& raw_surface) const
        {
            std::vector<VoxelGrid::GRID_INDEX> static_surface;
            // This may be larger than the actual surface we'll extrac
            static_surface.reserve(raw_surface.size());
            for (auto itr = raw_surface.begin(); itr != raw_surface.end(); ++itr)
            {
                const VoxelGrid::GRID_INDEX& index = itr->first;
                const uint8_t value = itr->second;
                if (value == 1)
                {
                    static_surface.push_back(index);
                }
            }
            // Try to reclaim the unnecessary vector capacity
            static_surface.shrink_to_fit();
            return static_surface;
        }

        std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t> ConvertToDynamicSurface(const std::vector<VoxelGrid::GRID_INDEX>& static_surface) const
        {
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t> dynamic_surface(static_surface.size());
            for (size_t idx = 0; idx < static_surface.size(); idx++)
            {
                const VoxelGrid::GRID_INDEX& grid_index = static_surface[idx];
                dynamic_surface[grid_index] = 1u;
            }
            return dynamic_surface;
        }

        std::unordered_map<VoxelGrid::GRID_INDEX, size_t> BuildSurfaceIndexMap(const std::vector<VoxelGrid::GRID_INDEX>& static_surface) const
        {
            std::unordered_map<VoxelGrid::GRID_INDEX, size_t> dynamic_surface(static_surface.size());
            for (size_t idx = 0; idx < static_surface.size(); idx++)
            {
                const VoxelGrid::GRID_INDEX& current_index = static_surface[idx];
                dynamic_surface[current_index] = idx;
            }
            return dynamic_surface;
        }

        std::pair<int32_t, int32_t> ComputeHolesInSurface(const uint32_t component, const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface, const bool verbose) const
        {
            // We have a list of all voxels with an exposed surface face
            // We loop through this list of voxels, and convert each voxel
            // into 8 vertices (the corners), which we individually check:
            //
            // First - we check to see if the vertex has already been
            // evaluated
            //
            // Second - we check if the vertex is actually on the surface
            // (make sure at least one of the three adjacent vertices is
            // exposed)
            //
            // Third - insert into hashtable of surface vertices
            //
            // Once we have completed this process, we loop back through
            // the hashtable of surface vertices and compute the number
            // of distance-1 neighboring surface vertices (we do this by
            // checking each of the six potential neighbor vertices) and
            // keep a running count of all vertices with 3, 5, and 6
            // neighbors.
            //
            // Once we have evaluated all the neighbors of all surface
            // vertices, we count the number of holes in the grid using
            // the formula from Chen and Rong, "Linear Time Recognition
            // Algorithms for Topological Invariants in 3D":
            //
            // #holes = 1 + (M5 + 2 * M6 - M3) / 8
            //
            // where M5 is the number of vertices with 5 neighbors,
            // M6 is the number of vertices with 6 neighbors, and
            // M3 is the number of vertices with 3 neighbors
            //
            // Storage for surface vertices
            // Compute a hint for initial surface vertex hashmap size
            // expected # of surface vertices
            // surface cells * 8
        #ifdef ENABLE_UNORDERED_MAP_SIZE_HINTS
            size_t surface_vertices_size_hint = surface.size() * 8;
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t> surface_vertices(surface_vertices_size_hint);
        #else
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t> surface_vertices;
        #endif
            // Loop through all the surface voxels and extract surface vertices
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>::const_iterator surface_itr;
            for (surface_itr = surface.begin(); surface_itr != surface.end(); ++surface_itr)
            {
                const VoxelGrid::GRID_INDEX& current_index = surface_itr->first;
                // First, grab all six neighbors from the grid
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyzm1 = GetImmutable(current_index.x, current_index.y, current_index.z - 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyzp1 = GetImmutable(current_index.x, current_index.y, current_index.z + 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xym1z = GetImmutable(current_index.x, current_index.y - 1, current_index.z);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xyp1z = GetImmutable(current_index.x, current_index.y + 1, current_index.z);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1yz = GetImmutable(current_index.x - 1, current_index.y, current_index.z);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1yz = GetImmutable(current_index.x + 1, current_index.y, current_index.z);
                // Generate all 8 vertices for the current voxel, check if an adjacent vertex is on the surface, and insert it if so
                // First, check the (-,-,-) vertex
                if (component != xyzm1.first.component || component != xym1z.first.component || component != xm1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex1(current_index.x, current_index.y, current_index.z);
                    surface_vertices[vertex1] = 1;
                }
                // Second, check the (-,-,+) vertex
                if (component != xyzp1.first.component || component != xym1z.first.component || component != xm1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex2(current_index.x, current_index.y, current_index.z + 1);
                    surface_vertices[vertex2] = 1;
                }
                // Third, check the (-,+,-) vertex
                if (component != xyzm1.first.component || component != xyp1z.first.component || component != xm1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex3(current_index.x, current_index.y + 1, current_index.z);
                    surface_vertices[vertex3] = 1;
                }
                // Fourth, check the (-,+,+) vertex
                if (component != xyzp1.first.component || component != xyp1z.first.component || component != xm1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex4(current_index.x, current_index.y + 1, current_index.z + 1);
                    surface_vertices[vertex4] = 1;
                }
                // Fifth, check the (+,-,-) vertex
                if (component != xyzm1.first.component || component != xym1z.first.component || component != xp1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex5(current_index.x + 1, current_index.y, current_index.z);
                    surface_vertices[vertex5] = 1;
                }
                // Sixth, check the (+,-,+) vertex
                if (component != xyzp1.first.component || component != xym1z.first.component || component != xp1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex6(current_index.x + 1, current_index.y, current_index.z + 1);
                    surface_vertices[vertex6] = 1;
                }
                // Seventh, check the (+,+,-) vertex
                if (component != xyzm1.first.component || component != xyp1z.first.component || component != xp1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex7(current_index.x + 1, current_index.y + 1, current_index.z);
                    surface_vertices[vertex7] = 1;
                }
                // Eighth, check the (+,+,+) vertex
                if (component != xyzp1.first.component || component != xyp1z.first.component || component != xp1yz.first.component)
                {
                    VoxelGrid::GRID_INDEX vertex8(current_index.x + 1, current_index.y + 1, current_index.z + 1);
                    surface_vertices[vertex8] = 1;
                }
            }
            if (verbose)
            {
                std::cerr << "Surface with " << surface.size() << " voxels has " << surface_vertices.size() << " surface vertices" << std::endl;
            }
            // Iterate through the surface vertices and count the neighbors of each vertex
            int32_t M3 = 0;
            int32_t M5 = 0;
            int32_t M6 = 0;
            // Store the connectivity of each vertex
            // Compute a hint for initial vertex connectivity hashmap size
            // real # of surface vertices
            // surface vertices
            size_t vertex_connectivity_size_hint = surface_vertices.size();
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t> vertex_connectivity(vertex_connectivity_size_hint);
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>::iterator surface_vertices_itr;
            for (surface_vertices_itr = surface_vertices.begin(); surface_vertices_itr != surface_vertices.end(); ++surface_vertices_itr)
            {
                VoxelGrid::GRID_INDEX key = surface_vertices_itr->first;
                VoxelGrid::GRID_INDEX value = key;
                // Insert into the connectivity map
                vertex_connectivity[key] = 0b00000000;
                // Check the six edges from the current vertex and count the number of exposed edges
                // (an edge is exposed if the at least one of the four surrounding voxels is not part
                // of the current component)
                int32_t edge_count = 0;
                // First, get the 8 voxels that surround the current vertex
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1ym1zm1 = GetImmutable(value.x - 1, value.y - 1, value.z - 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1ym1zp1 = GetImmutable(value.x - 1, value.y - 1, value.z + 0);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1yp1zm1 = GetImmutable(value.x - 1, value.y + 0, value.z - 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xm1yp1zp1 = GetImmutable(value.x - 1, value.y + 0, value.z + 0);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1ym1zm1 = GetImmutable(value.x + 0, value.y - 1, value.z - 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1ym1zp1 = GetImmutable(value.x + 0, value.y - 1, value.z + 0);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1yp1zm1 = GetImmutable(value.x + 0, value.y + 0, value.z - 1);
                std::pair<const TAGGED_OBJECT_COLLISION_CELL&, bool> xp1yp1zp1 = GetImmutable(value.x + 0, value.y + 0, value.z + 0);
                // Check the "z- down" edge
                if (component != xm1ym1zm1.first.component || component != xm1yp1zm1.first.component || component != xp1ym1zm1.first.component || component != xp1yp1zm1.first.component)
                {
                    if (!(component != xm1ym1zm1.first.component && component != xm1yp1zm1.first.component && component != xp1ym1zm1.first.component && component != xp1yp1zm1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00000001;
                    }
                }
                // Check the "z+ up" edge
                if (component != xm1ym1zp1.first.component || component != xm1yp1zp1.first.component || component != xp1ym1zp1.first.component || component != xp1yp1zp1.first.component)
                {
                    if (!(component != xm1ym1zp1.first.component && component != xm1yp1zp1.first.component && component != xp1ym1zp1.first.component && component != xp1yp1zp1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00000010;
                    }
                }
                // Check the "y- right" edge
                if (component != xm1ym1zm1.first.component || component != xm1ym1zp1.first.component || component != xp1ym1zm1.first.component || component != xp1ym1zp1.first.component)
                {
                    if (!(component != xm1ym1zm1.first.component && component != xm1ym1zp1.first.component && component != xp1ym1zm1.first.component && component != xp1ym1zp1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00000100;
                    }
                }
                // Check the "y+ left" edge
                if (component != xm1yp1zm1.first.component || component != xm1yp1zp1.first.component || component != xp1yp1zm1.first.component || component != xp1yp1zp1.first.component)
                {
                    if (!(component != xm1yp1zm1.first.component && component != xm1yp1zp1.first.component && component != xp1yp1zm1.first.component && component != xp1yp1zp1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00001000;
                    }
                }
                // Check the "x- back" edge
                if (component != xm1ym1zm1.first.component || component != xm1ym1zp1.first.component || component != xm1yp1zm1.first.component || component != xm1yp1zp1.first.component)
                {
                    if (!(component != xm1ym1zm1.first.component && component != xm1ym1zp1.first.component && component != xm1yp1zm1.first.component && component != xm1yp1zp1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00010000;
                    }
                }
                // Check the "x+ front" edge
                if (component != xp1ym1zm1.first.component || component != xp1ym1zp1.first.component || component != xp1yp1zm1.first.component || component != xp1yp1zp1.first.component)
                {
                    if (!(component != xp1ym1zm1.first.component && component != xp1ym1zp1.first.component && component != xp1yp1zm1.first.component && component != xp1yp1zp1.first.component))
                    {
                        edge_count++;
                        vertex_connectivity[key] |= 0b00100000;
                    }
                }
                // Increment M counts
                if (edge_count == 3)
                {
                    M3++;
                }
                else if (edge_count == 5)
                {
                    M5++;
                }
                else if (edge_count == 6)
                {
                    M6++;
                }
            }
            // Check to see if the set of vertices is connected. If not, our object contains void(s)
            int32_t number_of_surfaces = ComputeConnectivityOfSurfaceVertices(vertex_connectivity);
            int32_t number_of_voids = number_of_surfaces - 1;
            // Compute the number of holes in the surface
            int32_t raw_number_of_holes = 1 + ((M5 + (2 * M6) - M3) / 8);
            int32_t number_of_holes = raw_number_of_holes + number_of_voids;
            if (verbose)
            {
                std::cout << "Processing surface with M3 = " << M3 << " M5 = " << M5 << " M6 = " << M6 << " holes = " << number_of_holes << " surfaces = " << number_of_surfaces << " voids = " << number_of_voids << std::endl;
            }
            return std::pair<int32_t, int32_t>(number_of_holes, number_of_voids);
        }

        int32_t ComputeConnectivityOfSurfaceVertices(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface_vertex_connectivity) const
        {
            int32_t connected_components = 0;
            int64_t processed_vertices = 0;
            // Compute a hint for initial vertex components hashmap size
            // real # of surface vertices
            // surface vertices
            size_t vertex_components_size_hint = surface_vertex_connectivity.size();
            std::unordered_map<VoxelGrid::GRID_INDEX, int32_t> vertex_components(vertex_components_size_hint);
            // Iterate through the vertices
            std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>::const_iterator surface_vertices_itr;
            for (surface_vertices_itr = surface_vertex_connectivity.begin(); surface_vertices_itr != surface_vertex_connectivity.end(); ++surface_vertices_itr)
            {
                VoxelGrid::GRID_INDEX key = surface_vertices_itr->first;
                VoxelGrid::GRID_INDEX location = key;
                //const uint8_t& connectivity = surface_vertices_itr->second.second;
                // First, check if the vertex has already been marked
                if (vertex_components[key] > 0)
                {
                    continue;
                }
                else
                {
                    // If not, we start marking a new connected component
                    connected_components++;
                    // Make the working queue
                    std::list<VoxelGrid::GRID_INDEX> working_queue;
                    // Make a hash table to store queued indices (so we don't repeat work)
                    // Compute a hint for initial queued hashtable hashmap size
                    // If we assume that most object surfaces are, in fact, intact, then the first (and only)
                    // queued_hashtable will need to store an entry for every vertex on the surface.
                    // real # of surface vertices
                    // surface vertices
                    size_t queued_hashtable_size_hint = surface_vertex_connectivity.size();
                    std::unordered_map<VoxelGrid::GRID_INDEX, int8_t> queued_hashtable(queued_hashtable_size_hint);
                    // Add the current point
                    working_queue.push_back(location);
                    queued_hashtable[key] = 1;
                    // Keep track of the number of vertices we've processed
                    int64_t component_processed_vertices = 0;
                    // Loop from the queue
                    while (working_queue.size() > 0)
                    {
                        // Get the top of the working queue
                        VoxelGrid::GRID_INDEX current_vertex = working_queue.front();
                        working_queue.pop_front();
                        component_processed_vertices++;
                        vertex_components[current_vertex] = connected_components;
                        // Check the six possibly-connected vertices and add them to the queue if they are connected
                        // Get the connectivity of our index
                        uint8_t connectivity = surface_vertex_connectivity.at(current_vertex);
                        // Go through the neighbors
                        if ((connectivity & 0b00000001) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x, current_vertex.y, current_vertex.z - 1);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                        if ((connectivity & 0b00000010) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x, current_vertex.y, current_vertex.z + 1);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                        if ((connectivity & 0b00000100) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x, current_vertex.y - 1, current_vertex.z);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                        if ((connectivity & 0b00001000) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x, current_vertex.y + 1, current_vertex.z);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                        if ((connectivity & 0b00010000) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x - 1, current_vertex.y, current_vertex.z);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                        if ((connectivity & 0b00100000) > 0)
                        {
                            // Try to add the vertex
                            VoxelGrid::GRID_INDEX connected_vertex(current_vertex.x + 1, current_vertex.y, current_vertex.z);
                            // We only add if we haven't already processed it
                            if (queued_hashtable[connected_vertex] <= 0)
                            {
                                queued_hashtable[connected_vertex] = 1;
                                working_queue.push_back(connected_vertex);
                            }
                        }
                    }
                    processed_vertices += component_processed_vertices;
                    if (processed_vertices == (int64_t)surface_vertex_connectivity.size())
                    {
                        break;
                    }
                }
            }
            return connected_components;
        }

        inline std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> ExtractSignedDistanceField(const float oob_value, const std::vector<uint32_t>& objects_to_use) const
        {
            // To make this faster, we put the objects to use into a map
            std::map<uint32_t, uint8_t> object_use_map;
            for (size_t idx = 0; idx < objects_to_use.size(); idx++)
            {
                object_use_map[objects_to_use[idx]] = 1;
            }
            // Make the SDF
            SignedDistanceField new_sdf(GetOriginTransform(), frame_, GetResolution(), GetXSize(), GetYSize(), GetZSize(), oob_value);
            std::vector<VoxelGrid::GRID_INDEX> filled;
            std::vector<VoxelGrid::GRID_INDEX> free;
            for (int64_t x_index = 0; x_index < new_sdf.GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < new_sdf.GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < new_sdf.GetNumZCells(); z_index++)
                    {
                        VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
                        const TAGGED_OBJECT_COLLISION_CELL& stored = GetImmutable(x_index, y_index, z_index).first;
                        // If it matches an object to use OR there are no objects supplied
                        if ((object_use_map[stored.object_id] == 1) || (objects_to_use.size() == 0))
                        {
                            if (stored.occupancy > 0.5)
                            {
                                // Mark as filled
                                filled.push_back(current_index);
                            }
                            else
                            {
                                // Mark as free space
                                free.push_back(current_index);
                            }
                        }
                        else
                        {
                            // Mark as free space
                            free.push_back(current_index);
                        }
                    }
                }
            }
            // Make two distance fields (one for distance to filled voxels, one for distance to free voxels
            DistanceField filled_distance_field = BuildDistanceField(filled);
            DistanceField free_distance_field = BuildDistanceField(free);
            // Generate the SDF
            double max_distance = -INFINITY;
            double min_distance = INFINITY;
            for (int64_t x_index = 0; x_index < filled_distance_field.GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < filled_distance_field.GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < filled_distance_field.GetNumZCells(); z_index++)
                    {
                        double distance1 = sqrt(filled_distance_field.GetImmutable(x_index, y_index, z_index).first.distance_square) * new_sdf.GetResolution();
                        double distance2 = sqrt(free_distance_field.GetImmutable(x_index, y_index, z_index).first.distance_square) * new_sdf.GetResolution();
                        double distance = distance1 - distance2;
                        if (distance > max_distance)
                        {
                            max_distance = distance;
                        }
                        if (distance < min_distance)
                        {
                            min_distance = distance;
                        }
                        new_sdf.Set(x_index, y_index, z_index, distance);
                    }
                }
            }
            std::pair<double, double> extrema(max_distance, min_distance);
            return std::pair<SignedDistanceField, std::pair<double, double>>(new_sdf, extrema);
        }

        VoxelGrid::VoxelGrid<std::vector<uint32_t>> ComputeConvexRegions(const double max_check_radius) const;

        void GrowConvexRegion(const VoxelGrid::GRID_INDEX& start_index, VoxelGrid::VoxelGrid<std::vector<uint32_t>>& region_grid, const double max_check_radius, const uint32_t current_convex_region) const;

        EigenHelpers::VectorVector3d GenerateRayPrimitiveVectors(const uint32_t number_of_rays, const double cone_angle) const
        {
            (void)(number_of_rays);
            (void)(cone_angle);
            // Uniformly between [0.0, 2PI)
            std::vector<double> phis = {0.0, M_PI_4, M_PI_2, M_PI_4 * 3.0, M_PI, M_PI_4 * 5.0, M_PI_2 * 3.0, M_PI_4 * 7.0};
            // Uniformly between [cos(cone angle), 1.0]
            std::vector<double> zs = {0.55, 0.65, 0.75, 0.85, 0.95};
            // Sample the rays inside the cone
            EigenHelpers::VectorVector3d ray_primitive_vectors;
            for (size_t phidx = 0; phidx < phis.size(); phidx++)
            {
                for (size_t zdx = 0; zdx < zs.size(); zdx++)
                {
                    const double phi = phis[phidx];
                    const double z = zs[zdx];
                    // Compute the vector
                    const double x = sqrt(1.0 - (z * z)) * cos(phi);
                    const double y = sqrt(1.0 - (z * z)) * sin(phi);
                    Eigen::Vector3d raw_cone_vector(x, y, z);
                    ray_primitive_vectors.push_back(raw_cone_vector / raw_cone_vector.norm());
                }
            }
            ray_primitive_vectors.push_back(Eigen::Vector3d::UnitZ());
            return ray_primitive_vectors;
        }

        std::pair<std::vector<size_t>, std::vector<size_t>> CastSingleRay(const std::unordered_map<VoxelGrid::GRID_INDEX, size_t>& surface_index_map, const VoxelGrid::GRID_INDEX& current_surface_index, const Eigen::Vector3d& ray_unit_vector) const
        {
            std::map<size_t, uint8_t> line_of_sight_indices;
            std::map<size_t, uint8_t> non_line_of_sight_indices;
            Eigen::Vector3d current_location = EigenHelpers::StdVectorDoubleToEigenVector3d(GridIndexToLocation(current_surface_index.x, current_surface_index.y, current_surface_index.z));
            const double ray_step_size = GetResolution() * 0.5;
            const Eigen::Vector3d ray_step_vector = ray_unit_vector * ray_step_size;
            bool in_grid = true;
            bool collided = false;
            // Step along the ray vector until we run off the side of the grid
            while (in_grid)
            {
                // Grab the index corresponding to our location
                std::vector<int64_t> raw_current_index = LocationToGridIndex(current_location.x(), current_location.y(), current_location.z());
                if (raw_current_index.size() == 3)
                {
                    VoxelGrid::GRID_INDEX current_index(raw_current_index[0], raw_current_index[1], raw_current_index[2]);
                    //std::cout << "CVGX " << PrettyPrint::PrettyPrint(current_index) << std::endl;
                    // Are we in the surface?
                    auto found_itr = surface_index_map.find(current_index);
                    if (found_itr != surface_index_map.end())
                    {
                        //std::cout << "+++++ On surface +++++" << std::endl;
                        const VoxelGrid::GRID_INDEX& found_index = found_itr->first;
                        const size_t found_surface_index = found_itr->second;
                        // If we are not the current surface index
                        if (!(found_index == current_surface_index))
                        {
                            if ((line_of_sight_indices[found_surface_index] == 1) || (collided == false))
                            {
                                line_of_sight_indices[found_surface_index] = 1u;
                                collided = true;
                            }
                            else
                            {
                                non_line_of_sight_indices[found_surface_index] = 1u;
                            }
                        }
                    }
                    // Either way, we take a step
                    current_location += ray_step_vector;
                }
                // We're not in the grid any more
                else
                {
                    in_grid = false;
                }
            }
            //std::cout << "LOS [map] " << PrettyPrint::PrettyPrint(line_of_sight_indices) << std::endl;
            //std::cout << "NLOS [map] " << PrettyPrint::PrettyPrint(non_line_of_sight_indices) << std::endl;
            // Get the vectors of indices
            std::vector<size_t> line_of_sight_indices_vector;
            line_of_sight_indices_vector.reserve(line_of_sight_indices.size());
            for (auto itr = line_of_sight_indices.begin(); itr != line_of_sight_indices.end(); ++itr)
            {
                if (itr->second == 1u)
                {
                    line_of_sight_indices_vector.push_back(itr->first);
                }
            }
            std::vector<size_t> non_line_of_sight_indices_vector;
            non_line_of_sight_indices_vector.reserve(non_line_of_sight_indices.size());
            for (auto itr = non_line_of_sight_indices.begin(); itr != non_line_of_sight_indices.end(); ++itr)
            {
                if (itr->second == 1u)
                {
                    non_line_of_sight_indices_vector.push_back(itr->first);
                }
            }
            // Shrink to fit as needed
            line_of_sight_indices_vector.shrink_to_fit();
            non_line_of_sight_indices_vector.shrink_to_fit();
            //std::cout << "LOS [vec] " << PrettyPrint::PrettyPrint(line_of_sight_indices_vector) << std::endl;
            //std::cout << "NLOS [vec] " << PrettyPrint::PrettyPrint(non_line_of_sight_indices_vector) << std::endl;
            return std::pair<std::vector<size_t>, std::vector<size_t>>(line_of_sight_indices_vector, non_line_of_sight_indices_vector);
        }

        std::pair<std::vector<size_t>, std::vector<size_t>> PerformRayCasting(const sdf_tools::SignedDistanceField& sdf, const std::unordered_map<VoxelGrid::GRID_INDEX, size_t>& surface_index_map, const VoxelGrid::GRID_INDEX& current_surface_index, const EigenHelpers::VectorVector3d& ray_primitive_vectors) const
        {
            std::vector<size_t> line_of_sight_indices;
            std::vector<size_t> non_line_of_sight_indices;
            // Are we inside an object?
            bool inside_object = false;
            const float occupancy = GetImmutable(current_surface_index).first.occupancy;
            if (occupancy < 0.5f)
            {
                inside_object = false;
            }
            else if (occupancy > 0.5f)
            {
                inside_object = true;
            }
            else
            {
                // LOL NOPE
                assert(occupancy != 0.5f);
            }
            // Get the current gradient
            std::vector<double> raw_gradient = sdf.GetGradient(current_surface_index.x, current_surface_index.y, current_surface_index.z, true);
            Eigen::Vector3d raw_gradient_vector = EigenHelpers::StdVectorDoubleToEigenVector3d(raw_gradient);
            // Turn the gradient into an inverse surface normal vector (i.e. points inwards)
            assert(raw_gradient_vector.norm() > 0.0);
            // In free space, the vector already points in the correct direction
            Eigen::Vector3d surface_normal = raw_gradient_vector / raw_gradient_vector.norm();
            // If we're inside an object, flip the vector
            if (inside_object)
            {
                surface_normal = surface_normal * -1.0;
            }
            // Compute the pointing rotation using the cross product
            Eigen::Vector3d base_normalized_vector = Eigen::Vector3d::UnitZ();
            Eigen::Vector3d cross_product = base_normalized_vector.cross(surface_normal);
            double dot_product = base_normalized_vector.dot(surface_normal);
            double angle = acos(dot_product);
            // Make the rotation
            Eigen::AngleAxisd surface_normal_rotation(angle, cross_product);
            Eigen::Matrix3d surface_normal_rotation_matrix(surface_normal_rotation);
            // Safety check
            Eigen::Vector3d check_vector = surface_normal_rotation_matrix * base_normalized_vector;
            check_vector = check_vector / check_vector.norm();
            //std::cout << "Surface normal " << PrettyPrint::PrettyPrint(surface_normal) << std::endl;
            //std::cout << "Check vector " << PrettyPrint::PrettyPrint(check_vector) << std::endl;
            //assert(EigenHelpers::CloseEnough(check_vector, surface_normal, 0.01));
            // Perform raycasting for each ray
            for (size_t idx = 0; idx < ray_primitive_vectors.size(); idx++)
            {
                // Get the real ray unit vector
                const Eigen::Vector3d& ray_primitive_vector = ray_primitive_vectors[idx];
                const Eigen::Vector3d ray_unit_vector = surface_normal_rotation_matrix * ray_primitive_vector;
                //std::cout << "Ray vector " << PrettyPrint::PrettyPrint(ray_unit_vector) << std::endl;
                // Cast the single ray
                std::pair<std::vector<size_t>, std::vector<size_t>> single_raycast_results = CastSingleRay(surface_index_map, current_surface_index, ray_unit_vector);
                //std::cout << "Raycasting results - " << PrettyPrint::PrettyPrint(single_raycast_results) << std::endl;
                // Store the results
                line_of_sight_indices.insert(line_of_sight_indices.end(), single_raycast_results.first.begin(), single_raycast_results.first.end());
                non_line_of_sight_indices.insert(non_line_of_sight_indices.end(), single_raycast_results.second.begin(), single_raycast_results.second.end());
            }
            return std::pair<std::vector<size_t>, std::vector<size_t>>(line_of_sight_indices, non_line_of_sight_indices);
        }

        //std::pair<Eigen::MatrixXd, std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>> ComputeSparseLineOfSight(const std::vector<VoxelGrid::GRID_INDEX>& static_surface, const uint32_t number_of_rays, const double cone_angle) const
        Eigen::MatrixXd ComputeSparseLineOfSight(const std::vector<VoxelGrid::GRID_INDEX>& static_surface, const uint32_t number_of_rays, const double cone_angle) const
        {
            // Make our signed distance field
            sdf_tools::SignedDistanceField sdf = ExtractSignedDistanceField(INFINITY, std::vector<uint32_t>()).first;
            // Make a dynamic surface map
            std::unordered_map<VoxelGrid::GRID_INDEX, size_t> dynamic_surface_map = BuildSurfaceIndexMap(static_surface);
            // Make the ray primitive vectors
            EigenHelpers::VectorVector3d ray_primitive_vectors = GenerateRayPrimitiveVectors(number_of_rays, cone_angle);
            // Make containers
//            std::vector<Eigen::Triplet<double>> line_of_sight;
//            std::vector<Eigen::Triplet<double>> line_of_sight_not;
            Eigen::MatrixXd line_of_sight_matrix = Eigen::MatrixXd::Zero(static_surface.size(), static_surface.size());
            // Go through the surface
            for (size_t idx = 0; idx < static_surface.size(); idx++)
            {
                const VoxelGrid::GRID_INDEX& current_surface_index = static_surface[idx];
                // Perform raycasting
                // Returns a vector of line-of-sight surface indices, and a vector of non-line-of-sight surface indices
                const std::pair<std::vector<size_t>, std::vector<size_t>> raycast_results = PerformRayCasting(sdf, dynamic_surface_map, current_surface_index, ray_primitive_vectors);
                // Update the matrices
                for (size_t sdx = 0; sdx < raycast_results.first.size(); sdx++)
                {
                    const size_t surface_index = raycast_results.first[sdx];
//                    line_of_sight.push_back(Eigen::Triplet<double>(idx, surface_index, 1.0));
//                    line_of_sight.push_back(Eigen::Triplet<double>(surface_index, idx, 1.0));
                    line_of_sight_matrix(idx, surface_index) = 1.0;
                    line_of_sight_matrix(surface_index, idx) = 1.0;
                }
                for (size_t sdx = 0; sdx < raycast_results.second.size(); sdx++)
                {
                    const size_t surface_index = raycast_results.second[sdx];
//                    line_of_sight_not.push_back(Eigen::Triplet<double>(idx, surface_index, 1.0));
//                    line_of_sight_not.push_back(Eigen::Triplet<double>(surface_index, idx, 1.0));
                    line_of_sight_matrix(idx, surface_index) = 0.0;
                    line_of_sight_matrix(surface_index, idx) = 0.0;
                }
                // Add ourselves to the line of sight
//                line_of_sight.push_back(Eigen::Triplet<double>(idx, idx, 1.0));
                line_of_sight_matrix(idx, idx) = 1.0;
            }
//            // Put together the sparse matrices
//            Eigen::SparseMatrix<double> line_of_sight_sparse;
//            line_of_sight_sparse.setFromTriplets(line_of_sight.begin(), line_of_sight.end());
//            Eigen::SparseMatrix<double> line_of_sight_not_sparse;
//            line_of_sight_not_sparse.setFromTriplets(line_of_sight_not.begin(), line_of_sight_not.end());
            // Return them all
            //return std::pair<Eigen::MatrixXd, std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>>(line_of_sight_matrix, std::pair<Eigen::SparseMatrix<double>, Eigen::SparseMatrix<double>>(line_of_sight_sparse, line_of_sight_not_sparse));
            return line_of_sight_matrix;
        }

        std::pair<Eigen::VectorXd, Eigen::MatrixXd> ExtractKLargestEigenvaluesAndEigenvectors(const Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType& raw_eigenvalues, const Eigen::EigenSolver<Eigen::MatrixXd>::EigenvectorsType& raw_eigenvectors, const uint32_t num_values) const
        {
            assert((int64_t)num_values <= raw_eigenvalues.size());
            // Collect the eigenvalue + index pairs
            std::vector<std::pair<double, size_t>> eigenvalue_index_pairs(raw_eigenvalues.size());
            for (size_t idx = 0; idx < eigenvalue_index_pairs.size(); idx++)
            {
                const double current_eigenvalue = raw_eigenvalues((long)idx).real();
                eigenvalue_index_pairs[idx].first = current_eigenvalue;
                eigenvalue_index_pairs[idx].second = idx;
            }
            // Sort the eigenvalue/index pairs by eigenvalue
            std::function<bool(const std::pair<double, size_t>&, const std::pair<double, size_t>&)> compare_fn = [] (const std::pair<double, size_t>& p1, const std::pair<double, size_t>& p2) { if (p1.first < p2.first) { return true; } else { return false; } };
            // Sorts them in ascending order
            std::sort(eigenvalue_index_pairs.begin(), eigenvalue_index_pairs.end(), compare_fn);
            // Now, we extract the last num_values eigenvalues and eigenvectors and put them into a vector + matrix
            // Each column of this matrix is an eigenvector, so the # of rows corresponds to the number of datapoints, and the # of columns is the number of clusters
            Eigen::VectorXd k_largest_eigenvalues = Eigen::VectorXd::Zero(num_values);
            Eigen::MatrixXd k_largest_eigenvectors = Eigen::MatrixXd::Zero(raw_eigenvalues.rows(), num_values);
            for (uint32_t num_value = 0; num_value < num_values; num_value++)
            {
                // Get the index of the ascending-sorted vector
                const int64_t eigenvalue_index = (eigenvalue_index_pairs.size() - 1) - num_value;
                // Get the index corresponding to the num_valueth-largest eigenvalue
                const double current_eigenvalue = eigenvalue_index_pairs[eigenvalue_index].first;
                k_largest_eigenvalues(num_value) = current_eigenvalue;
                const size_t eigenvector_index = eigenvalue_index_pairs[eigenvalue_index].second;
                // Grab the corresponding eigenvector
                const Eigen::VectorXcd current_eigenvector = raw_eigenvectors.col((int64_t)eigenvector_index);
                // Copy it over into the real world
                Eigen::VectorXd real_current_eigenvector = Eigen::VectorXd::Zero(current_eigenvector.size());
                for (int64_t vdx = 0; vdx < real_current_eigenvector.size(); vdx++)
                {
                    real_current_eigenvector(vdx) = current_eigenvector(vdx).real();
                }
                // Set it in the matrix
                k_largest_eigenvectors.col(num_value) = real_current_eigenvector;
            }
            // Return
            return std::pair<Eigen::VectorXd, Eigen::MatrixXd>(k_largest_eigenvalues, k_largest_eigenvectors);
        }

        std::vector<uint32_t> PerformKMeansSpectralClustering(const Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType& raw_eigenvalues, const Eigen::EigenSolver<Eigen::MatrixXd>::EigenvectorsType& raw_eigenvectors, const uint32_t num_clusters) const
        {
            assert(num_clusters > 0);
            // Grab the num_clusters largest eigenvalues & corresponding eigenvectors
            // Each column of this matrix is an eigenvector, so the # of rows corresponds to the number of datapoints, and the # of columns is the number of clusters
            Eigen::MatrixXd k_largest_eigenvectors = ExtractKLargestEigenvaluesAndEigenvectors(raw_eigenvalues, raw_eigenvectors, num_clusters).second;
            // Convert them into datapoints for kmeans clustering
            std::vector<Eigen::VectorXd> clustering_data(k_largest_eigenvectors.rows());
            for (size_t datapoint_index = 0; datapoint_index < clustering_data.size(); datapoint_index++)
            {
                Eigen::VectorXd datapoint = Eigen::VectorXd::Zero(num_clusters);
                for (int64_t datavalue_index = 0; datavalue_index < (int64_t)num_clusters; datavalue_index++)
                {
                    const double raw_value = k_largest_eigenvectors(datapoint_index, datavalue_index);
                    datapoint(datavalue_index) = raw_value;
                }
                clustering_data[datapoint_index] = datapoint;
            }
            std::function<double(const Eigen::VectorXd&, const Eigen::VectorXd&)> distance_fn = [] (const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) { return EigenHelpers::Distance(v1, v2); };
            std::function<Eigen::VectorXd(const std::vector<Eigen::VectorXd>&)> average_fn = [] (const std::vector<Eigen::VectorXd>& data) { return EigenHelpers::AverageEigenVectorXd(data); };
            std::vector<uint32_t> cluster_labels = simple_kmeans_clustering::SimpleKMeansClustering::Cluster(clustering_data, distance_fn, average_fn, num_clusters, true);
            return cluster_labels;
        }

        double ComputeConvexityMetric(const Eigen::MatrixXd& los_matrix, const std::vector<uint32_t>& cluster_labels) const
        {
            const double alpha = 1.0; // This was used in the paper
            assert(los_matrix.rows() == los_matrix.cols());
            assert(cluster_labels.size() == (size_t)los_matrix.rows());
            const uint32_t num_clusters = *std::max_element(cluster_labels.begin(), cluster_labels.end()) + 1;
            double total_convexity_metric = 0.0;
            for (uint32_t cluster = 0; cluster < num_clusters; cluster++)
            {
                double intravisible = 0.0;
                double interoccluded = 0.0;
                for (size_t idx = 0; idx < cluster_labels.size(); idx++)
                {
                    // We only care about elements in our own cluster
                    if (cluster_labels[idx] == cluster)
                    {
                        // Loop through our row in the LOS matrix
                        for (size_t other_idx = 0; other_idx < cluster_labels.size(); other_idx++)
                        {
                            // If the other element is part of our cluster AND is visible
                            if ((cluster_labels[other_idx] == cluster) && (los_matrix(idx, other_idx) == 1.0))
                            {
                                intravisible += 1.0; // This double counts, but we need to double count for convexity = 1 in the single cluster case
                            }
                            // If the other element is not part of our cluster AND is not visible
                            else if ((cluster_labels[other_idx] != cluster) && (los_matrix(idx, other_idx) == 0.0))
                            {
                                interoccluded += 1.0;
                            }
                        }
                    }
                }
                const double cluster_convexity_metric = intravisible + (alpha * interoccluded);
                total_convexity_metric += cluster_convexity_metric;
            }
            total_convexity_metric = total_convexity_metric / (double)(cluster_labels.size() * cluster_labels.size());
            return total_convexity_metric;
        }

        std::vector<std::vector<size_t>> ClusterSurfaceFromLOSMatrix(const Eigen::MatrixXd& los_matrix, const uint32_t max_num_clusters) const
        {
            assert(los_matrix.rows() == los_matrix.cols());
            assert(max_num_clusters > 0);
            // Perform spectral clustering
            // Compute degree matrix
            Eigen::MatrixXd degree_matrix_diagonals = los_matrix.rowwise().sum();
            Eigen::MatrixXd degree_matrix = Eigen::MatrixXd::Zero(los_matrix.rows(), los_matrix.cols());
            degree_matrix.diagonal() = degree_matrix_diagonals;
            // Compute the unnormalized Laplacian
            Eigen::MatrixXd unnormalized_laplacian = degree_matrix - los_matrix;
            // Compute Eigenvalues & Eigenvectors
            Eigen::EigenSolver<Eigen::MatrixXd> solver(unnormalized_laplacian);
            Eigen::EigenSolver<Eigen::MatrixXd>::EigenvalueType eigen_values = solver.eigenvalues();
            Eigen::EigenSolver<Eigen::MatrixXd>::EigenvectorsType eigen_vectors = solver.eigenvectors();
            // Perform k-means clustering over a range of # of clusters
            double best_convexity = 0.0;
            std::vector<uint32_t> best_clustering;
            for (uint32_t num_clusters = 1; num_clusters <= max_num_clusters; num_clusters++)
            {
                if (num_clusters > (uint32_t)los_matrix.rows())
                {
                    std::cerr << "Number of clusters is larger than elements in surface, stopping clustering process" << std::endl;
                    break;
                }
                std::vector<uint32_t> cluster_labels = PerformKMeansSpectralClustering(eigen_values, eigen_vectors, num_clusters);
                double convexity = ComputeConvexityMetric(los_matrix, cluster_labels);
                std::cout << "K-means clustering at " << num_clusters << " clusters with convexity " << convexity << std::endl;
                if (convexity > best_convexity)
                {
                    best_convexity = convexity;
                    best_clustering = cluster_labels;
                }
            }
            // Safety check
            assert(best_clustering.size() == (size_t)los_matrix.rows());
            // Turn the clustering labels into separate clusters
            const uint32_t num_clusters = *std::max_element(best_clustering.begin(), best_clustering.end()) + 1;
            std::vector<std::vector<size_t>> clustered_surfaces(num_clusters);
            for (size_t idx = 0; idx < best_clustering.size(); idx++)
            {
                const uint32_t cluster_label = best_clustering[idx];
                clustered_surfaces[cluster_label].push_back(idx);
            }
            return clustered_surfaces;
        }

        std::vector<std::vector<VoxelGrid::GRID_INDEX>> ComputeWeaklyConvexSurfaceSegments(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface, const uint32_t max_num_clusters) const
        {
            std::vector<VoxelGrid::GRID_INDEX> static_surface = ExtractStaticSurface(surface);
            Eigen::MatrixXd LOS_matrix = ComputeSparseLineOfSight(static_surface, 30u, (M_PI_2 * 0.5)); // For now, these values are actually ignored
            std::cout << "LOS matrix:\n" << LOS_matrix << std::endl;
            std::vector<std::vector<size_t>> convex_surface_segment_indices = ClusterSurfaceFromLOSMatrix(LOS_matrix, max_num_clusters);
            // Convert the 0-n indices into grid indices
            std::vector<std::vector<VoxelGrid::GRID_INDEX>> convex_surface_segments(convex_surface_segment_indices.size());
            for (size_t segment_idx = 0; segment_idx < convex_surface_segment_indices.size(); segment_idx++)
            {
                const std::vector<size_t>& current_segment_indices = convex_surface_segment_indices[segment_idx];
                convex_surface_segments[segment_idx].reserve(current_segment_indices.size());
                for (size_t index_idx = 0; index_idx < current_segment_indices.size(); index_idx++)
                {
                    const size_t current_segment_index = current_segment_indices[index_idx];
                    const VoxelGrid::GRID_INDEX& current_segment_surface_index = static_surface[current_segment_index];
                    convex_surface_segments[segment_idx].push_back(current_segment_surface_index);
                }
                assert(convex_surface_segments[segment_idx].size() == current_segment_indices.size());
            }
            assert(convex_surface_segment_indices.size() == convex_surface_segments.size());
            return convex_surface_segments;
        }

        std::map<uint32_t, uint32_t> UpdateConvexSegments()
        {
            // Some day, we will do real work here. Until then, this is a dummy that does nothing
            convex_segments_valid_ = true;
            // Return a map of object_id to # of convex segments in the object
            std::map<uint32_t, uint32_t> convex_segment_counts;
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        const TAGGED_OBJECT_COLLISION_CELL& cell = GetImmutable(x_index, y_index, z_index).first;
                        const uint32_t cell_object_id = cell.object_id;
                        const std::vector<uint32_t> cell_convex_segments = cell.GetListOfConvexSegments();
                        if (cell_convex_segments.size() > 0)
                        {
                            const uint32_t max_segment_number = *std::max_element(cell_convex_segments.begin(), cell_convex_segments.end());
                            if (max_segment_number > convex_segment_counts[cell_object_id])
                            {
                                convex_segment_counts[cell_object_id] = max_segment_number;
                            }
                        }
                    }
                }
            }
            return convex_segment_counts;
        }

        std::map<uint32_t, sdf_tools::SignedDistanceField> MakeObjectSDFs(const std::vector<uint32_t>& object_ids) const
        {
            std::map<uint32_t, sdf_tools::SignedDistanceField> per_object_sdfs;
            for (size_t idx = 0; idx < object_ids.size(); idx++)
            {
                const uint32_t object_id = object_ids[idx];
                per_object_sdfs[object_id] = ExtractSignedDistanceField(std::numeric_limits<double>::infinity(), std::vector<uint32_t>{object_id}).first;
            }
            return per_object_sdfs;
        }

        std::map<uint32_t, sdf_tools::SignedDistanceField> MakeObjectSDFs() const
        {
            std::map<uint32_t, uint32_t> object_id_map;
            for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
                    {
                        const TAGGED_OBJECT_COLLISION_CELL& cell = GetImmutable(x_index, y_index, z_index).first;
                        const uint32_t cell_object_id = cell.object_id;
                        if (cell_object_id > 0)
                        {
                            object_id_map[cell_object_id] = 1u;
                        }
                    }
                }
            }
            return MakeObjectSDFs(arc_helpers::GetKeys(object_id_map));
        }

        visualization_msgs::Marker ExportForDisplay(const float alpha, const std::vector<uint32_t>& objects_to_draw=std::vector<uint32_t>()) const;

        visualization_msgs::Marker ExportForDisplay(const std::map<uint32_t, std_msgs::ColorRGBA>& object_color_map=std::map<uint32_t, std_msgs::ColorRGBA>()) const;

        visualization_msgs::Marker ExportContourOnlyForDisplay(const float alpha, const std::vector<uint32_t>& objects_to_draw=std::vector<uint32_t>()) const;

        visualization_msgs::Marker ExportContourOnlyForDisplay(const std::map<uint32_t, std_msgs::ColorRGBA>& object_color_map=std::map<uint32_t, std_msgs::ColorRGBA>()) const;

        visualization_msgs::Marker ExportForDisplayOccupancyOnly(const std_msgs::ColorRGBA& collision_color, const std_msgs::ColorRGBA& free_color, const std_msgs::ColorRGBA& unknown_color) const;

        visualization_msgs::Marker ExportConnectedComponentsForDisplay(bool color_unknown_components) const;

        visualization_msgs::Marker ExportConvexSegmentForDisplay(const uint32_t object_id, const uint32_t convex_segment) const;

        visualization_msgs::Marker ExportSurfaceForDisplay(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface, const std_msgs::ColorRGBA& surface_color) const;
    };
}

#endif // TAGGED_OBJECT_COLLISION_MAP_HPP
