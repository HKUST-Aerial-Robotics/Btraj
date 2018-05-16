#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <Eigen/Geometry>
#include <visualization_msgs/Marker.h>
#include <arc_utilities/eigen_helpers.hpp>
#include <arc_utilities/voxel_grid.hpp>
#include <sdf_tools/SDF.h>

#ifndef SDF_HPP
#define SDF_HPP

inline std::vector<uint8_t> FloatToBinary(float value)
{
    uint32_t binary_value = 0;
    memcpy(&binary_value, &value, sizeof(uint32_t));
    std::vector<uint8_t> binary(4);
    // Copy byte 1, least-significant byte
    binary[3] = binary_value & 0x000000ff;
    // Copy byte 2
    binary_value = binary_value >> 8;
    binary[2] = binary_value & 0x000000ff;
    // Copy byte 3
    binary_value = binary_value >> 8;
    binary[1] = binary_value & 0x000000ff;
    // Copy byte 4, most-significant byte
    binary_value = binary_value >> 8;
    binary[0] = binary_value & 0x000000ff;
    return binary;
}

inline float FloatFromBinary(std::vector<uint8_t>& binary)
{
    if (binary.size() != 4)
    {
        std::cerr << "Binary value is not 4 bytes" << std::endl;
        return NAN;
    }
    else
    {
        uint32_t binary_value = 0;
        // Copy in byte 4, most-significant byte
        binary_value = binary_value | binary[0];
        binary_value = binary_value << 8;
        // Copy in byte 3
        binary_value = binary_value | binary[1];
        binary_value = binary_value << 8;
        // Copy in byte 2
        binary_value = binary_value | binary[2];
        binary_value = binary_value << 8;
        // Copy in byte 1, least-significant byte
        binary_value = binary_value | binary[3];
        // Convert binary to float and store
        float field_value = 0.0;
        memcpy(&field_value, &binary_value, sizeof(float));
        return field_value;
    }
}

namespace sdf_tools
{
    class SignedDistanceField
    {
    protected:

        VoxelGrid::VoxelGrid<float> distance_field_;
        std::string frame_;
        bool initialized_;
        bool locked_;

        std::vector<uint8_t> GetInternalBinaryRepresentation(const std::vector<float> &field_data);

        std::vector<float> UnpackFieldFromBinaryRepresentation(std::vector<uint8_t>& binary);

        /*
         * You *MUST* provide valid indices to this function, hence why it is protected (there are safe wrappers available - use them!)
         */
        void FollowGradientsToLocalMaximaUnsafe(VoxelGrid::VoxelGrid<Eigen::Vector3d>& watershed_map, const int64_t x_index, const int64_t y_index, const int64_t z_index) const;

    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

        inline SignedDistanceField(std::string frame, double resolution, double x_size, double y_size, double z_size, float OOB_value) : initialized_(true), locked_(false)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<float> new_field(resolution, x_size, y_size, z_size, OOB_value);
            distance_field_ = new_field;
        }

        inline SignedDistanceField(Eigen::Affine3d origin_transform, std::string frame, double resolution, double x_size, double y_size, double z_size, float OOB_value) : initialized_(true), locked_(false)
        {
            frame_ = frame;
            VoxelGrid::VoxelGrid<float> new_field(origin_transform, resolution, x_size, y_size, z_size, OOB_value);
            distance_field_ = new_field;
        }

        inline SignedDistanceField() : initialized_(false), locked_(false) {}

        inline bool IsInitialized() const
        {
            return initialized_;
        }

        inline bool IsLocked() const
        {
            return locked_;
        }

        inline void Lock()
        {
            locked_ = true;
        }

        inline void Unlock()
        {
            locked_ = false;
        }

        inline float Get(const double x, const double y, const double z) const
        {
            return distance_field_.GetImmutable(x, y, z).first;
        }

        inline float Get3d(const Eigen::Vector3d& location) const
        {
            return distance_field_.GetImmutable3d(location).first;
        }

        inline float Get4d(const Eigen::Vector4d& location) const
        {
            return distance_field_.GetImmutable4d(location).first;
        }

        inline float Get(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return distance_field_.GetImmutable(x_index, y_index, z_index).first;
        }

        inline std::pair<float, bool> GetSafe(const double x, const double y, const double z) const
        {
            return distance_field_.GetImmutable(x, y, z);
        }

        inline std::pair<float, bool> GetSafe3d(const Eigen::Vector3d& location) const
        {
            return distance_field_.GetImmutable3d(location);
        }

        inline std::pair<float, bool> GetSafe4d(const Eigen::Vector4d& location) const
        {
            return distance_field_.GetImmutable4d(location);
        }

        inline std::pair<float, bool> GetSafe(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return distance_field_.GetImmutable(x_index, y_index, z_index);
        }

        /*
         * Setter functions MUST be used carefully - If you arbitrarily change SDF values, it is not a proper SDF any more!
         *
         * Use of these functions can be prevented by calling SignedDistanceField::Lock() on the SDF, at which point these functions
         * will fail with a warning printed to std_err.
         */
        inline bool Set(const double x, const double y, const double z, float value)
        {
            if (!locked_)
            {
                return distance_field_.SetValue(x, y, z, value);
            }
            else
            {
                std::cerr << "Attempt to set value in locked SDF" << std::endl;
                return false;
            }
        }

        inline bool Set3d(const Eigen::Vector3d& location, float value)
        {
            if (!locked_)
            {
                return distance_field_.SetValue3d(location, value);
            }
            else
            {
                std::cerr << "Attempt to set value in locked SDF" << std::endl;
                return false;
            }
        }

        inline bool Set4d(const Eigen::Vector4d& location, float value)
        {
            if (!locked_)
            {
                return distance_field_.SetValue4d(location, value);
            }
            else
            {
                std::cerr << "Attempt to set value in locked SDF" << std::endl;
                return false;
            }
        }

        inline bool Set(const int64_t x_index, const int64_t y_index, const int64_t z_index, const float value)
        {
            if (!locked_)
            {
                return distance_field_.SetValue(x_index, y_index, z_index, value);
            }
            else
            {
                std::cerr << "Attempt to set value in locked SDF" << std::endl;
                return false;
            }
        }

        inline bool Set(const VoxelGrid::GRID_INDEX& index, const float value)
        {
            if (!locked_)
            {
                return distance_field_.SetValue(index, value);
            }
            else
            {
                std::cerr << "Attempt to set value in locked SDF" << std::endl;
                return false;
            }
        }

        inline bool CheckInBounds3d(const Eigen::Vector3d& location) const
        {
            return distance_field_.GetImmutable3d(location).second;
        }

        inline bool CheckInBounds4d(const Eigen::Vector4d& location) const
        {
            return distance_field_.GetImmutable4d(location).second;
        }

        inline bool CheckInBounds(const double x, const double y, const double z) const
        {
            return distance_field_.GetImmutable(x, y, z).second;
        }

        inline bool CheckInBounds(const VoxelGrid::GRID_INDEX& index) const
        {
            return distance_field_.GetImmutable(index.x, index.y, index.z).second;
        }

        inline bool CheckInBounds(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return distance_field_.GetImmutable(x_index, y_index, z_index).second;
        }

        inline double GetXSize() const
        {
            return distance_field_.GetXSize();
        }

        inline double GetYSize() const
        {
            return distance_field_.GetYSize();
        }

        inline double GetZSize() const
        {
            return distance_field_.GetZSize();
        }

        inline double GetResolution() const
        {
            return distance_field_.GetCellSizes()[0];
        }

        inline float GetOOBValue() const
        {
            return distance_field_.GetDefaultValue();
        }

        inline int64_t GetNumXCells() const
        {
            return distance_field_.GetNumXCells();
        }

        inline int64_t GetNumYCells() const
        {
            return distance_field_.GetNumYCells();
        }

        inline int64_t GetNumZCells() const
        {
            return distance_field_.GetNumZCells();
        }

    protected:

        inline std::pair<Eigen::Vector3d, double> GetPrimaryComponentsVector(const Eigen::Vector3d& raw_vector) const
        {
            if (std::abs(raw_vector.x()) > std::abs(raw_vector.y()) && std::abs(raw_vector.x()) > std::abs(raw_vector.z()))
            {
                if (raw_vector.x() >= 0.0)
                {
                    return std::make_pair(Eigen::Vector3d(GetResolution() * 0.5, 0.0, 0.0), GetResolution() * 0.5);
                }
                else
                {
                    return std::make_pair(Eigen::Vector3d(GetResolution() * -0.5, 0.0, 0.0), GetResolution() * 0.5);
                }
            }
            else if (std::abs(raw_vector.y()) > std::abs(raw_vector.x()) && std::abs(raw_vector.y()) > std::abs(raw_vector.z()))
            {
                if (raw_vector.y() >= 0.0)
                {
                    return std::make_pair(Eigen::Vector3d(0.0, GetResolution() * 0.5, 0.0), GetResolution() * 0.5);
                }
                else
                {
                    return std::make_pair(Eigen::Vector3d(0.0, GetResolution() * -0.5, 0.0), GetResolution() * 0.5);
                }
            }
            else if (std::abs(raw_vector.z()) > std::abs(raw_vector.x()) && std::abs(raw_vector.z()) > std::abs(raw_vector.y()))
            {
                if (raw_vector.z() >= 0.0)
                {
                    return std::make_pair(Eigen::Vector3d(0.0, 0.0, GetResolution() * 0.5), GetResolution() * 0.5);
                }
                else
                {
                    return std::make_pair(Eigen::Vector3d(0.0, 0.0, GetResolution() * -0.5), GetResolution() * 0.5);
                }
            }
            else if (std::abs(raw_vector.x()) == std::abs(raw_vector.y()))
            {
                const Eigen::Vector3d temp_vector(raw_vector.x(), raw_vector.y(), 0.0);
                return std::make_pair((temp_vector / (temp_vector.norm())) * std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0), std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0));
            }
            else if (std::abs(raw_vector.y()) == std::abs(raw_vector.z()))
            {
                const Eigen::Vector3d temp_vector(0.0, raw_vector.y(), raw_vector.x());
                return std::make_pair((temp_vector / (temp_vector.norm())) * std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0), std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0));
            }
            else if (std::abs(raw_vector.x()) == std::abs(raw_vector.z()))
            {
                const Eigen::Vector3d temp_vector(raw_vector.x(), 0.0, raw_vector.z());
                return std::make_pair((temp_vector / (temp_vector.norm())) * std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0), std::sqrt((GetResolution() * GetResolution() * 0.25) * 2.0));
            }
            else
            {
                return std::make_pair((raw_vector / (raw_vector.norm())) * std::sqrt((GetResolution() * GetResolution() * 0.25) * 3.0), std::sqrt((GetResolution() * GetResolution() * 0.25) * 3.0));
            }
        }

        inline double ComputeAxisMatch(const double axis_value, const double check_value) const
        {
            if ((axis_value >= 0.0) == (check_value >= 0.0))
            {
                return std::abs(check_value - axis_value);
            }
            else
            {
                return -std::abs(check_value - axis_value);
            }
        }

        inline Eigen::Vector3d GetBestMatchSurfaceVector(const Eigen::Vector3d& possible_surfaces_vector, const Eigen::Vector3d& center_to_location_vector) const
        {
            const Eigen::Vector3d location_rejected_on_possible = EigenHelpers::VectorRejection(possible_surfaces_vector, center_to_location_vector);
            // Find the axis with the best-match components
            const double x_axis_match = ComputeAxisMatch(possible_surfaces_vector.x(), location_rejected_on_possible.x());
            const double y_axis_match = ComputeAxisMatch(possible_surfaces_vector.y(), location_rejected_on_possible.y());
            const double z_axis_match = ComputeAxisMatch(possible_surfaces_vector.z(), location_rejected_on_possible.z());
            if ((x_axis_match > y_axis_match) && (x_axis_match > z_axis_match))
            {
                return Eigen::Vector3d(possible_surfaces_vector.x(), 0.0, 0.0);
            }
            else if ((y_axis_match > x_axis_match) && (y_axis_match > z_axis_match))
            {
                return Eigen::Vector3d(0.0, possible_surfaces_vector.y(), 0.0);
            }
            else if ((z_axis_match > x_axis_match) && (z_axis_match > y_axis_match))
            {
                return Eigen::Vector3d(0.0, 0.0, possible_surfaces_vector.z());
            }
            else
            {
                assert(false);
                return possible_surfaces_vector;
            }
        }

        /**
         * @brief GetPrimaryEntrySurfaceVector Estimates the real distance of the provided point, comparing it with the cell center location and gradient vector
         * @param boundary_direction_vector
         * @param center_to_location_vector
         * @return vector from center of voxel to primary entry surface, and magnitude of that vector
         */
        inline std::pair<Eigen::Vector3d, double> GetPrimaryEntrySurfaceVector(const Eigen::Vector3d& boundary_direction_vector, const Eigen::Vector3d& center_to_location_vector) const
        {
            if (boundary_direction_vector.squaredNorm() > std::numeric_limits<double>::epsilon())
            {
                const std::pair<Eigen::Vector3d, double> primary_components_vector_query = GetPrimaryComponentsVector(boundary_direction_vector);
                // If the cell is on a surface
                if (primary_components_vector_query.second == (GetResolution() * 0.5))
                {
                    return primary_components_vector_query;
                }
                // If the cell is on an edge or surface
                else
                {
                    // Pick the best-match of the two/three exposed surfaces
                    return std::make_pair(GetBestMatchSurfaceVector(primary_components_vector_query.first, center_to_location_vector), GetResolution() * 0.5);
                }
            }
            else
            {
                return GetPrimaryComponentsVector(center_to_location_vector);
            }
        }

        inline double EstimateDistanceInternal(const double x, const double y, const double z, const int64_t x_idx, const int64_t y_idx, const int64_t z_idx) const
        {
            const std::vector<double> cell_center = GridIndexToLocation(x_idx, y_idx, z_idx);
            const Eigen::Vector3d cell_center_to_location_vector(x - cell_center[0], y - cell_center[1], z - cell_center[2]);
            const double nominal_sdf_distance = (double)distance_field_.GetImmutable(x_idx, y_idx, z_idx).first;

            // Determine vector from "entry surface" to center of voxel
            // TODO: Needs special handling if there's no gradient to work with
            const std::vector<double> raw_gradient = GetGradient(x_idx, y_idx, z_idx, true);
            const Eigen::Vector3d gradient = EigenHelpers::StdVectorDoubleToEigenVector3d(raw_gradient);
            const Eigen::Vector3d direction_to_boundary = (nominal_sdf_distance >= 0.0) ? -gradient : gradient;
            const std::pair<Eigen::Vector3d, double> entry_surface_information = GetPrimaryEntrySurfaceVector(direction_to_boundary, cell_center_to_location_vector);
            const Eigen::Vector3d& entry_surface_vector = entry_surface_information.first;
            const double minimum_distance_magnitude = entry_surface_information.second;

            // Adjust for calculating distance to boundary of voxels instead of center of voxels
            const double center_adjusted_nominal_distance = (nominal_sdf_distance >= 0.0) ? nominal_sdf_distance - (GetResolution() * 0.5) : nominal_sdf_distance + (GetResolution() * 0.5);
            const double minimum_adjusted_distance = arc_helpers::SpreadValue(center_adjusted_nominal_distance, -minimum_distance_magnitude, 0.0, minimum_distance_magnitude);

            // Account for target location being not at the exact center of the voxel
            const double raw_distance_adjustment = EigenHelpers::VectorProjection(entry_surface_vector, cell_center_to_location_vector).norm();
            const double real_distance_adjustment = (minimum_adjusted_distance >= 0.0) ? -raw_distance_adjustment: raw_distance_adjustment;
            const double final_adjusted_distance = minimum_adjusted_distance + real_distance_adjustment;

            // Perform minimum distance thresholding and error checking
            // TODO: do we need to address this magic number somehow?
            if (std::abs(final_adjusted_distance) < GetResolution() * 0.001)
            {
                return 0.0;
            }
            if ((minimum_adjusted_distance >= 0.0) == (final_adjusted_distance >= 0.0))
            {
                return final_adjusted_distance;
            }
            else
            {
                std::cerr << "Center adjusted nominal distance " << minimum_adjusted_distance << " final adjusted_distance " << final_adjusted_distance << std::endl;
                assert(false && "Mismatched minimum and final adjusted distance signs");
            }
        }

    public:

        inline std::pair<double, bool> EstimateDistance(const double x, const double y, const double z) const
        {
            return EstimateDistance4d(Eigen::Vector4d(x, y, z, 1.0));
        }

        inline std::pair<double, bool> EstimateDistance3d(const Eigen::Vector3d& location) const
        {
            const std::vector<int64_t> indices = LocationToGridIndex3d(location);
            if (indices.size() == 3)
            {
                return std::make_pair(EstimateDistanceInternal(location.x(), location.y(), location.z(), indices[0], indices[1], indices[2]), true);
            }
            else
            {
                return std::make_pair((double)distance_field_.GetOOBValue(), false);
            }
        }

        inline std::pair<double, bool> EstimateDistance4d(const Eigen::Vector4d& location) const
        {
            const std::vector<int64_t> indices = LocationToGridIndex4d(location);
            if (indices.size() == 3)
            {
                return std::make_pair(EstimateDistanceInternal(location(0), location(1), location(2), indices[0], indices[1], indices[2]), true);
            }
            else
            {
                return std::make_pair((double)distance_field_.GetOOBValue(), false);
            }
        }

        inline std::vector<double> GetGradient(const double x, const double y, const double z, const bool enable_edge_gradients=false) const
        {
            return GetGradient4d(Eigen::Vector4d(x, y, z, 1.0), enable_edge_gradients);
        }

        inline std::vector<double> GetGradient3d(const Eigen::Vector3d& location, const bool enable_edge_gradients=false) const
        {
            const std::vector<int64_t> indices = LocationToGridIndex3d(location);
            if (indices.size() == 3)
            {
                return GetGradient(indices[0], indices[1], indices[2], enable_edge_gradients);
            }
            else
            {
                return std::vector<double>();
            }
        }

        inline std::vector<double> GetGradient4d(const Eigen::Vector4d& location, const bool enable_edge_gradients=false) const
        {
            const std::vector<int64_t> indices = LocationToGridIndex4d(location);
            if (indices.size() == 3)
            {
                return GetGradient(indices[0], indices[1], indices[2], enable_edge_gradients);
            }
            else
            {
                return std::vector<double>();
            }
        }

        inline std::vector<double> GetGradient(const VoxelGrid::GRID_INDEX& index, const bool enable_edge_gradients=false) const
        {
            return GetGradient(index.x, index.y, index.z, enable_edge_gradients);
        }

        inline std::vector<double> GetGradient(const int64_t x_index, const int64_t y_index, const int64_t z_index, const bool enable_edge_gradients=false) const
        {
            // Make sure the index is inside bounds
            if ((x_index >= 0) && (y_index >= 0) && (z_index >= 0) && (x_index < GetNumXCells()) && (y_index < GetNumYCells()) && (z_index < GetNumZCells()))
            {
                // Make sure the index we're trying to query is one cell in from the edge
                if ((x_index > 0) && (y_index > 0) && (z_index > 0) && (x_index < (GetNumXCells() - 1)) && (y_index < (GetNumYCells() - 1)) && (z_index < (GetNumZCells() - 1)))
                {
                    double inv_twice_resolution = 1.0 / (2.0 * GetResolution());
                    double gx = (Get(x_index + 1, y_index, z_index) - Get(x_index - 1, y_index, z_index)) * inv_twice_resolution;
                    double gy = (Get(x_index, y_index + 1, z_index) - Get(x_index, y_index - 1, z_index)) * inv_twice_resolution;
                    double gz = (Get(x_index, y_index, z_index + 1) - Get(x_index, y_index, z_index - 1)) * inv_twice_resolution;
                    return std::vector<double>{gx, gy, gz};
                }
                // If we're on the edge, handle it specially
                else if (enable_edge_gradients)
                {
                    // Get the "best" indices we can use
                    int64_t low_x_index = std::max((int64_t)0, x_index - 1);
                    int64_t high_x_index = std::min(GetNumXCells() - 1, x_index + 1);
                    int64_t low_y_index = std::max((int64_t)0, y_index - 1);
                    int64_t high_y_index = std::min(GetNumYCells() - 1, y_index + 1);
                    int64_t low_z_index = std::max((int64_t)0, z_index - 1);
                    int64_t high_z_index = std::min(GetNumZCells() - 1, z_index + 1);
                    // Compute the axis increments
                    double x_increment = (high_x_index - low_x_index) * GetResolution();
                    double y_increment = (high_y_index - low_y_index) * GetResolution();
                    double z_increment = (high_z_index - low_z_index) * GetResolution();
                    // Compute the gradients for each axis - by default these are zero
                    double gx = 0.0;
                    double gy = 0.0;
                    double gz = 0.0;
                    // Only if the increments are non-zero do we compute the gradient of an axis
                    if (x_increment > 0.0)
                    {
                        double inv_x_increment = 1.0 / x_increment;
                        double high_x_value = Get(high_x_index, y_index, z_index);
                        double low_x_value = Get(low_x_index, y_index, z_index);
                        // Compute the gradient
                        gx = (high_x_value - low_x_value) * inv_x_increment;
                    }
                    if (y_increment > 0.0)
                    {
                        double inv_y_increment = 1.0 / y_increment;
                        double high_y_value = Get(x_index, high_y_index, z_index);
                        double low_y_value = Get(x_index, low_y_index, z_index);
                        // Compute the gradient
                        gy = (high_y_value - low_y_value) * inv_y_increment;
                    }
                    if (z_increment > 0.0)
                    {
                        double inv_z_increment = 1.0 / z_increment;
                        double high_z_value = Get(x_index, y_index, high_z_index);
                        double low_z_value = Get(x_index, y_index, low_z_index);
                        // Compute the gradient
                        gz = (high_z_value - low_z_value) * inv_z_increment;
                    }
                    // Assemble and return the computed gradient
                    return std::vector<double>{gx, gy, gz};
                }
                // Edge gradients disabled, return no gradient
                else
                {
                    return std::vector<double>();
                }
            }
            // If we're out of bounds, return no gradient
            else
            {
                return std::vector<double>();
            }
        }

        inline Eigen::Vector3d ProjectOutOfCollision(const double x, const double y, const double z, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            const Eigen::Vector4d result = ProjectOutOfCollision4d(Eigen::Vector4d(x, y, z, 1.0), stepsize_multiplier);
            return result.head<3>();
        }

        inline Eigen::Vector3d ProjectOutOfCollisionToMinimumDistance(const double x, const double y, const double z, const double minimum_distance, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            const Eigen::Vector4d result = ProjectOutOfCollisionToMinimumDistance4d(Eigen::Vector4d(x, y, z, 1.0), minimum_distance, stepsize_multiplier);
            return result.head<3>();
        }

        inline Eigen::Vector3d ProjectOutOfCollision3d(const Eigen::Vector3d& location, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            return ProjectOutOfCollision(location.x(), location.y(), location.z(), stepsize_multiplier);
        }

        inline Eigen::Vector3d ProjectOutOfCollisionToMinimumDistance3d(const Eigen::Vector3d& location, const double minimum_distance, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            return ProjectOutOfCollisionToMinimumDistance(location.x(), location.y(), location.z(), minimum_distance, stepsize_multiplier);
        }

        inline Eigen::Vector4d ProjectOutOfCollision4d(const Eigen::Vector4d& location, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            return ProjectOutOfCollisionToMinimumDistance4d(location, 0.0, stepsize_multiplier);
        }

        inline Eigen::Vector4d ProjectOutOfCollisionToMinimumDistance4d(const Eigen::Vector4d& location, const double minimum_distance, const double stepsize_multiplier = 1.0 / 10.0) const
        {
            Eigen::Vector4d mutable_location = location;
            const bool enable_edge_gradients = true;

            double sdf_dist = EstimateDistance4d(mutable_location).first;
            if (sdf_dist < minimum_distance && CheckInBounds4d(location))
            {
                while (sdf_dist < minimum_distance)
                {
                    const std::vector<double> gradient = GetGradient4d(mutable_location, enable_edge_gradients);
                    const Eigen::Vector3d grad_eigen = EigenHelpers::StdVectorDoubleToEigenVector3d(gradient);

                    assert(grad_eigen.norm() > GetResolution() / 4.0); // Sanity check
                    mutable_location.head<3>() += grad_eigen.normalized() * GetResolution() * stepsize_multiplier;

                    sdf_dist = EstimateDistance4d(mutable_location).first;
                }
            }

            return mutable_location;
        }

        inline const Eigen::Affine3d& GetOriginTransform() const
        {
            return distance_field_.GetOriginTransform();
        }

        inline const Eigen::Affine3d& GetInverseOriginTransform() const
        {
            return distance_field_.GetInverseOriginTransform();
        }

        inline std::string GetFrame() const
        {
            return frame_;
        }

        inline std::vector<int64_t> LocationToGridIndex3d(const Eigen::Vector3d& location) const
        {
            return distance_field_.LocationToGridIndex3d(location);
        }

        inline std::vector<int64_t> LocationToGridIndex4d(const Eigen::Vector4d& location) const
        {
            return distance_field_.LocationToGridIndex4d(location);
        }

        inline std::vector<int64_t> LocationToGridIndex(const double x, const double y, const double z) const
        {
            return distance_field_.LocationToGridIndex(x, y, z);
        }

        inline std::vector<double> GridIndexToLocation(const VoxelGrid::GRID_INDEX& index) const
        {
            return distance_field_.GridIndexToLocation(index);
        }

        inline std::vector<double> GridIndexToLocation(const int64_t x_index, const int64_t y_index, const int64_t z_index) const
        {
            return distance_field_.GridIndexToLocation(x_index, y_index, z_index);
        }

        inline std::vector<double> GridIndexToLocation(std::vector<int64_t> index) const
        {
            return distance_field_.GridIndexToLocation(index[0], index[1], index[2]);
        }

        bool SaveToFile(const std::string& filepath);

        bool LoadFromFile(const std::string& filepath);

        sdf_tools::SDF GetMessageRepresentation();

        bool LoadFromMessageRepresentation(sdf_tools::SDF& message);

        visualization_msgs::Marker ExportForDisplay(float alpha = 0.01f) const;

        visualization_msgs::Marker ExportForDisplayCollisionOnly(float alpha = 0.01f) const;

        visualization_msgs::Marker ExportForDebug(float alpha = 0.5f) const;

        /*
         * The following function can be *VERY EXPENSIVE* to compute, since it performs gradient ascent across the SDF
         */
        VoxelGrid::VoxelGrid<Eigen::Vector3d> ComputeLocalMaximaMap() const;

        inline bool GradientIsEffectiveFlat(const Eigen::Vector3d& gradient) const
        {
            // A gradient is at a local maxima if the absolute value of all components (x,y,z) are less than 1/2 SDF resolution
            double half_resolution = GetResolution() * 0.5;
            if (fabs(gradient.x()) <= half_resolution && fabs(gradient.y()) <= half_resolution && fabs(gradient.z()) <= half_resolution)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        inline VoxelGrid::GRID_INDEX GetNextFromGradient(const VoxelGrid::GRID_INDEX& index, const Eigen::Vector3d& gradient) const
        {
            // Given the gradient, pick the "best fit" of the 26 neighboring points
            VoxelGrid::GRID_INDEX next_index = index;
            double half_resolution = GetResolution() * 0.5;
            if (gradient.x() > half_resolution)
            {
                next_index.x++;
            }
            else if (gradient.x() < -half_resolution)
            {
                next_index.x--;
            }
            if (gradient.y() > half_resolution)
            {
                next_index.y++;
            }
            else if (gradient.y() < -half_resolution)
            {
                next_index.y--;
            }
            if (gradient.z() > half_resolution)
            {
                next_index.z++;
            }
            else if (gradient.z() < -half_resolution)
            {
                next_index.z--;
            }
            return next_index;
        }
    };
}

#endif // SDF_HPP
