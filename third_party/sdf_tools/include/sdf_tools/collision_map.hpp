#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <Eigen/Geometry>
#include <visualization_msgs/Marker.h>
#include <arc_utilities/arc_helpers.hpp>
#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/sdf.hpp>
#include <sdf_tools/CollisionMap.h>

#include <eigen3/Eigen/Dense>

#ifndef COLLISION_MAP_HPP
#define COLLISION_MAP_HPP

#define ENABLE_UNORDERED_MAP_SIZE_HINTS

namespace sdf_tools
{
    struct COLLISION_CELL
    {
        float occupancy;
        uint32_t component;

        COLLISION_CELL() : occupancy(0.0), component(0) {}

        COLLISION_CELL(const float in_occupancy) : occupancy(in_occupancy), component(0) {}

        COLLISION_CELL(const float in_occupancy, const uint32_t in_component) : occupancy(in_occupancy), component(in_component) {}
    };

    inline std::vector<uint8_t> CollisionCellToBinary(const COLLISION_CELL& value)
    {
        std::vector<uint8_t> binary(sizeof(COLLISION_CELL));
        memcpy(&binary.front(), &value, sizeof(COLLISION_CELL));
        return binary;
    }

    inline COLLISION_CELL CollisionCellFromBinary(const std::vector<uint8_t>& binary)
    {
        if (binary.size() != sizeof(COLLISION_CELL))
        {
            std::cerr << "Binary value is not " << sizeof(COLLISION_CELL) << " bytes" << std::endl;
            return COLLISION_CELL(NAN, 0u);
        }
        else
        {
            COLLISION_CELL loaded;
            memcpy(&loaded, &binary.front(), sizeof(COLLISION_CELL));
            return loaded;
        }
    }

    class CollisionMapGrid
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
            uint32_t our_component = collision_field_.GetImmutable(x_index, y_index, z_index).first.component;
            // Check neighbor 1
            if (our_component != collision_field_.GetImmutable(x_index, y_index, z_index - 1).first.component)
            {
                return true;
            }
            // Check neighbor 2
            else if (our_component != collision_field_.GetImmutable(x_index, y_index, z_index + 1).first.component)
            {
                return true;
            }
            // Check neighbor 3
            else if (our_component != collision_field_.GetImmutable(x_index, y_index - 1, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 4
            else if (our_component != collision_field_.GetImmutable(x_index, y_index + 1, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 5
            else if (our_component != collision_field_.GetImmutable(x_index - 1, y_index, z_index).first.component)
            {
                return true;
            }
            // Check neighbor 6
            else if (our_component != collision_field_.GetImmutable(x_index + 1, y_index, z_index).first.component)
            {
                return true;
            }
            // If none of the faces are exposed, it's not a surface voxel
            return false;
        }

        typedef struct
        {
            uint32_t location[3];
            uint32_t closest_point[3];
            double distance_square;
            int32_t update_direction;
        } bucket_cell;

        typedef VoxelGrid::VoxelGrid<bucket_cell> DistanceField;

        inline DistanceField BuildDistanceField(const std::vector<VoxelGrid::GRID_INDEX>& points) const
        {
            // Make the DistanceField container
            bucket_cell default_cell;
            default_cell.distance_square = INFINITY;
            DistanceField distance_field(collision_field_.GetOriginTransform(), GetResolution(), collision_field_.GetXSize(), collision_field_.GetYSize(), collision_field_.GetZSize(), default_cell);
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

        VoxelGrid::VoxelGrid<COLLISION_CELL> collision_field_;
        uint32_t number_of_components_;
        std::string frame_;
        bool initialized_;
        bool components_valid_;

        std::vector<uint8_t> PackBinaryRepresentation(std::vector<COLLISION_CELL>& raw);

        std::vector<COLLISION_CELL> UnpackBinaryRepresentation(std::vector<uint8_t>& packed);

        int64_t MarkConnectedComponent(int64_t x_index, int64_t y_index, int64_t z_index, uint32_t connected_component);

    public:

        inline CollisionMapGrid(const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const COLLISION_CELL& default_value, const COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<COLLISION_CELL> new_field(resolution, x_size, y_size, z_size, default_value, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
        }

        inline CollisionMapGrid(const Eigen::Affine3d& origin_transform, const std::string& frame, const double resolution, const double x_size, double y_size, const double z_size, const COLLISION_CELL& default_value, const COLLISION_CELL& OOB_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<COLLISION_CELL> new_field(origin_transform, resolution, x_size, y_size, z_size, default_value, OOB_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
        }

        inline CollisionMapGrid(const std::string& frame, const double resolution, const double x_size, const double y_size, const double z_size, const COLLISION_CELL& OOB_default_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<COLLISION_CELL> new_field(resolution, x_size, y_size, z_size, OOB_default_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
        }

        inline CollisionMapGrid(const Eigen::Affine3d& origin_transform, const std::string& frame, const double resolution, const double x_size, double y_size, const double z_size, const COLLISION_CELL& OOB_default_value) : initialized_(true)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<COLLISION_CELL> new_field(origin_transform, resolution, x_size, y_size, z_size, OOB_default_value);
            collision_field_ = new_field;
            number_of_components_ = 0;
            components_valid_ = false;
        }

        inline CollisionMapGrid() : number_of_components_(0), initialized_(false), components_valid_(false) {}

        inline bool IsInitialized() const
        {
            return initialized_;
        }

        inline bool AreComponentsValid() const
        {
            return components_valid_;
        }

        inline std::pair<COLLISION_CELL, bool> Get3d(const Eigen::Vector3d& location) const
        {
            return collision_field_.GetImmutable3d(location);
        }

        inline std::pair<COLLISION_CELL, bool> Get4d(const Eigen::Vector4d& location) const
        {
            return collision_field_.GetImmutable4d(location);
        }

        inline std::pair<COLLISION_CELL, bool> Get(const double x, const double y, const double z) const
        {
            return collision_field_.GetImmutable(x, y, z);
        }

        inline std::pair<COLLISION_CELL, bool> Get(const VoxelGrid::GRID_INDEX& index) const
        {
            return collision_field_.GetImmutable(index);
        }

        inline std::pair<COLLISION_CELL, bool> Get(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return collision_field_.GetImmutable(x_index, y_index, z_index);
        }

        inline bool Set(const double x, const double y, const double z, COLLISION_CELL value)
        {
            components_valid_ = false;
            return collision_field_.SetValue(x, y, z, value);
        }

        void Set3d(const Eigen::Vector3d& location, COLLISION_CELL value)
        {
            //components_valid_ = false;
            collision_field_.SetValue3d(location, value);
        }

        inline bool Set4d(const Eigen::Vector4d& location, COLLISION_CELL value)
        {
            components_valid_ = false;
            return collision_field_.SetValue4d(location, value);
        }

        inline bool Set(const int64_t x_index, const int64_t y_index, const int64_t z_index, COLLISION_CELL value)
        {
            components_valid_ = false;
            return collision_field_.SetValue(x_index, y_index, z_index, value);
        }

        inline bool Set(const VoxelGrid::GRID_INDEX& index, COLLISION_CELL value)
        {
            components_valid_ = false;
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

        inline COLLISION_CELL GetDefaultValue() const
        {
            return collision_field_.GetDefaultValue();
        }

        inline COLLISION_CELL GetOOBValue() const
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

        inline std::vector<int64_t> LocationToGridIndex(double x, double y, double z) const
        {
            return collision_field_.LocationToGridIndex(x, y, z);
        }

        inline bool Inside(Eigen::Vector3i index) const
        {
            return collision_field_.IndexInBounds(index(0), index(1), index(2));
        }

        inline Eigen::Vector3i LocationToGridIndex(Eigen::Vector3d location) const
        {
            return collision_field_.LocationToGridIndex(location);
        }

        inline std::vector<double> GridIndexToLocation(int64_t x_index, int64_t y_index, int64_t z_index) const
        {
            return collision_field_.GridIndexToLocation(x_index, y_index, z_index);
        }

        inline Eigen::Vector3d GridIndexToLocation(Eigen::Vector3i index) const
        {
            return collision_field_.GridIndexToLocation(index);
        }

        bool SaveToFile(const std::string& filepath);

        bool LoadFromFile(const std::string &filepath);

        sdf_tools::CollisionMap GetMessageRepresentation();

        bool LoadFromMessageRepresentation(sdf_tools::CollisionMap& message);

        uint32_t UpdateConnectedComponents();

        std::map<uint32_t, std::pair<int32_t, int32_t>> ComputeComponentTopology(bool ignore_empty_components, bool recompute_connected_components, bool verbose);

        std::map<uint32_t, std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>> ExtractComponentSurfaces(const bool ignore_empty_components) const;

        std::pair<int32_t, int32_t> ComputeHolesInSurface(const uint32_t component, const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface, const bool verbose) const;

        int32_t ComputeConnectivityOfSurfaceVertices(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface_vertex_connectivity) const;

        inline std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> ExtractSignedDistanceField(const float oob_value) const
        {
            // Make the SDF
            SignedDistanceField new_sdf(collision_field_.GetOriginTransform(), frame_, GetResolution(), collision_field_.GetXSize(), collision_field_.GetYSize(), collision_field_.GetZSize(), oob_value);
            std::vector<VoxelGrid::GRID_INDEX> filled;
            std::vector<VoxelGrid::GRID_INDEX> free;
            for (int64_t x_index = 0; x_index < new_sdf.GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < new_sdf.GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < new_sdf.GetNumZCells(); z_index++)
                    {
                        VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
                        COLLISION_CELL stored = Get(x_index, y_index, z_index).first;
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
                        //double distance = distance1;
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

        inline DistanceField ExtractDistanceField(const float oob_value) const
        {
            // Make the SDF
            SignedDistanceField new_sdf(collision_field_.GetOriginTransform(), frame_, GetResolution(), collision_field_.GetXSize(), collision_field_.GetYSize(), collision_field_.GetZSize(), oob_value);
            std::vector<VoxelGrid::GRID_INDEX> filled;
            for (int64_t x_index = 0; x_index < new_sdf.GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < new_sdf.GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < new_sdf.GetNumZCells(); z_index++)
                    {
                        VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
                        if (Get(x_index, y_index, z_index).first.occupancy > 0.5)
                        {
                            // Mark as filled
                            filled.push_back(current_index);
                        }
                    }
                }
            }
            // Make two distance fields (one for distance to filled voxels, one for distance to free voxels
            DistanceField filled_distance_field = BuildDistanceField(filled);

            return filled_distance_field;
        }

        void RestMap()
        {
            // Reset components first
            for (int64_t x_index = 0; x_index < collision_field_.GetNumXCells(); x_index++)
            {
                for (int64_t y_index = 0; y_index < collision_field_.GetNumYCells(); y_index++)
                {
                    for (int64_t z_index = 0; z_index < collision_field_.GetNumZCells(); z_index++)
                    {
                        COLLISION_CELL free_cell(0.0);
                        collision_field_.SetValue(x_index, y_index, z_index, free_cell);
                    }
                }
            }
        }

        visualization_msgs::Marker ExportForDisplay(const std_msgs::ColorRGBA& collision_color, const std_msgs::ColorRGBA& free_color, const std_msgs::ColorRGBA& unknown_color) const;

        visualization_msgs::Marker ExportConnectedComponentsForDisplay(bool color_unknown_components) const;
    };
}

#endif // COLLISION_MAP_HPP
