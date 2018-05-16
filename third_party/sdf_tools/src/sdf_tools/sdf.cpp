#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <zlib.h>
#include <ros/ros.h>
#include <arc_utilities/eigen_helpers_conversions.hpp>
#include <arc_utilities/zlib_helpers.hpp>
#include <sdf_tools/sdf.hpp>
#include <sdf_tools/SDF.h>

using namespace sdf_tools;

std::vector<uint8_t> SignedDistanceField::GetInternalBinaryRepresentation(const std::vector<float>& field_data)
{
    std::vector<uint8_t> raw_binary_data(field_data.size() * 4);
    for (size_t field_index = 0, binary_index = 0; field_index < field_data.size(); field_index++, binary_index+=4)
    {
        // Convert the float at the current index into 4 bytes and store them
        float field_value = field_data[field_index];
        std::vector<uint8_t> binary_value = FloatToBinary(field_value);
        raw_binary_data[binary_index] = binary_value[0];
        raw_binary_data[binary_index + 1] = binary_value[1];
        raw_binary_data[binary_index + 2] = binary_value[2];
        raw_binary_data[binary_index + 3] = binary_value[3];
    }
    return raw_binary_data;
}

std::vector<float> SignedDistanceField::UnpackFieldFromBinaryRepresentation(std::vector<uint8_t>& binary)
{
    if ((binary.size() % 4) != 0)
    {
        std::cerr << "Invalid binary representation - length is not a multiple of 4" << std::endl;
        return std::vector<float>();
    }
    uint64_t data_size = binary.size() / 4;
    std::vector<float> field_data(data_size);
    for (size_t field_index = 0, binary_index = 0; field_index < field_data.size(); field_index++, binary_index+=4)
    {
        std::vector<uint8_t> binary_block{binary[binary_index], binary[binary_index + 1], binary[binary_index + 2], binary[binary_index + 3]};
        field_data[field_index] = FloatFromBinary(binary_block);
    }
    return field_data;
}

bool SignedDistanceField::SaveToFile(const std::string& filepath)
{
    // Convert to message representation
    sdf_tools::SDF message_rep = GetMessageRepresentation();
    // Save message to file
    try
    {
        std::ofstream output_file(filepath.c_str(), std::ios::out|std::ios::binary);
        uint32_t serialized_size = ros::serialization::serializationLength(message_rep);
        std::unique_ptr<uint8_t> ser_buffer(new uint8_t[serialized_size]);
        ros::serialization::OStream ser_stream(ser_buffer.get(), serialized_size);
        ros::serialization::serialize(ser_stream, message_rep);
        output_file.write((char*)ser_buffer.get(), serialized_size);
        output_file.close();
        return true;
    }
    catch (...)
    {
        return false;
    }
}

bool SignedDistanceField::LoadFromFile(const std::string &filepath)
{
    try
    {
        // Load message from file
        std::ifstream input_file(filepath.c_str(), std::ios::in|std::ios::binary);
        input_file.seekg(0, std::ios::end);
        std::streampos end = input_file.tellg();
        input_file.seekg(0, std::ios::beg);
        std::streampos begin = input_file.tellg();
        uint32_t serialized_size = end - begin;
        std::unique_ptr<uint8_t> deser_buffer(new uint8_t[serialized_size]);
        input_file.read((char*) deser_buffer.get(), serialized_size);
        ros::serialization::IStream deser_stream(deser_buffer.get(), serialized_size);
        sdf_tools::SDF new_message;
        ros::serialization::deserialize(deser_stream, new_message);
        // Load state from the message
        bool success = LoadFromMessageRepresentation(new_message);
        return success;
    }
    catch (...)
    {
        return false;
    }
}

sdf_tools::SDF SignedDistanceField::GetMessageRepresentation()
{
    sdf_tools::SDF message_rep;
    // Populate message
    message_rep.header.frame_id = frame_;
    const Eigen::Affine3d& origin_transform = distance_field_.GetOriginTransform();
    message_rep.origin_transform.translation.x = origin_transform.translation().x();
    message_rep.origin_transform.translation.y = origin_transform.translation().y();
    message_rep.origin_transform.translation.z = origin_transform.translation().z();
    const Eigen::Quaterniond origin_transform_rotation(origin_transform.rotation());
    message_rep.origin_transform.rotation.x = origin_transform_rotation.x();
    message_rep.origin_transform.rotation.y = origin_transform_rotation.y();
    message_rep.origin_transform.rotation.z = origin_transform_rotation.z();
    message_rep.origin_transform.rotation.w = origin_transform_rotation.w();
    message_rep.dimensions.x = distance_field_.GetXSize();
    message_rep.dimensions.y = distance_field_.GetYSize();
    message_rep.dimensions.z = distance_field_.GetZSize();
    message_rep.sdf_cell_size = GetResolution();
    message_rep.OOB_value = distance_field_.GetDefaultValue();
    message_rep.initialized = initialized_;
    message_rep.locked = locked_;
    const std::vector<float>& raw_data = distance_field_.GetRawData();
    std::vector<uint8_t> binary_data = GetInternalBinaryRepresentation(raw_data);
    message_rep.data = ZlibHelpers::CompressBytes(binary_data);
    return message_rep;
}

bool SignedDistanceField::LoadFromMessageRepresentation(sdf_tools::SDF& message)
{
    // Make a new voxel grid inside
    Eigen::Translation3d origin_translation(message.origin_transform.translation.x, message.origin_transform.translation.y, message.origin_transform.translation.z);
    Eigen::Quaterniond origin_rotation(message.origin_transform.rotation.w, message.origin_transform.rotation.x, message.origin_transform.rotation.y, message.origin_transform.rotation.z);
    Eigen::Affine3d origin_transform = origin_translation * origin_rotation;
    VoxelGrid::VoxelGrid<float> new_field(origin_transform, message.sdf_cell_size, message.dimensions.x, message.dimensions.y, message.dimensions.z, message.OOB_value);
    // Unpack the binary data
    std::vector<uint8_t> binary_data = ZlibHelpers::DecompressBytes(message.data);
    std::vector<float> unpacked = UnpackFieldFromBinaryRepresentation(binary_data);
    if (unpacked.empty())
    {
        std::cerr << "Unpack returned an empty SDF" << std::endl;
        return false;
    }
    bool success = new_field.SetRawData(unpacked);
    if (!success)
    {
        std::cerr << "Unable to set internal representation of the SDF" << std::endl;
        return false;
    }
    // Set it
    distance_field_ = new_field;
    frame_ = message.header.frame_id;
    initialized_ = message.initialized;
    locked_ = message.locked;
    return true;
}

visualization_msgs::Marker SignedDistanceField::ExportForDisplay(float alpha) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "sdf_display";
    display_rep.id = 1;
    display_rep.type = visualization_msgs::Marker::CUBE_LIST;
    display_rep.action = visualization_msgs::Marker::ADD;
    display_rep.lifetime = ros::Duration(0.0);
    display_rep.frame_locked = false;
    const Eigen::Affine3d base_transform = Eigen::Affine3d::Identity();
    display_rep.pose = EigenHelpersConversions::EigenAffine3dToGeometryPose(base_transform);
    display_rep.scale.x = GetResolution();
    display_rep.scale.y = GetResolution();
    display_rep.scale.z = GetResolution();
    // Add all the cells of the SDF to the message
    double min_distance = 0.0;
    double max_distance = 0.0;
    for (int64_t x_index = 0; x_index < distance_field_.GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < distance_field_.GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < distance_field_.GetNumZCells(); z_index++)
            {
                // Update minimum/maximum distance variables
                float distance = Get(x_index, y_index, z_index);
                if (distance < min_distance)
                {
                    min_distance = distance;
                }
                if (distance > max_distance)
                {
                    max_distance = distance;
                }
                // Convert SDF indices into a real-world location
                std::vector<double> location = distance_field_.GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                display_rep.points.push_back(new_point);
            }
        }
    }
    // Add colors for all the cells of the SDF to the message
    for (int64_t x_index = 0; x_index < distance_field_.GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < distance_field_.GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < distance_field_.GetNumZCells(); z_index++)
            {
                // Update minimum/maximum distance variables
                float distance = Get(x_index, y_index, z_index);
                std_msgs::ColorRGBA new_color;
                new_color.a = alpha;
                if (distance > 0.0)
                {
                    new_color.b = 0.0;
                    new_color.g = (fabs(distance / max_distance) * 0.8) + 0.2;
                    new_color.r = 0.0;
                }
                else if (distance < 0.0)
                {
                    new_color.b = 0.0;
                    new_color.g = 0.0;
                    new_color.r = (fabs(distance / min_distance) * 0.8) + 0.2;
                }
                else
                {
                    new_color.b = 1.0;
                    new_color.g = 0.0;
                    new_color.r = 0.0;
                }
                display_rep.colors.push_back(new_color);
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker SignedDistanceField::ExportForDisplayCollisionOnly(float alpha) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "sdf_display";
    display_rep.id = 1;
    display_rep.type = visualization_msgs::Marker::CUBE_LIST;
    display_rep.action = visualization_msgs::Marker::ADD;
    display_rep.lifetime = ros::Duration(0.0);
    display_rep.frame_locked = false;
    const Eigen::Affine3d base_transform = Eigen::Affine3d::Identity();
    display_rep.pose = EigenHelpersConversions::EigenAffine3dToGeometryPose(base_transform);
    display_rep.scale.x = GetResolution();
    display_rep.scale.y = GetResolution();
    display_rep.scale.z = GetResolution();
    // Add all the cells of the SDF to the message
    for (int64_t x_index = 0; x_index < distance_field_.GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < distance_field_.GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < distance_field_.GetNumZCells(); z_index++)
            {
                // Update minimum/maximum distance variables
                float distance = Get(x_index, y_index, z_index);
                if (distance <= 0.0)
                {
                    // Convert SDF indices into a real-world location
                    std::vector<double> location = distance_field_.GridIndexToLocation(x_index, y_index, z_index);
                    geometry_msgs::Point new_point;
                    new_point.x = location[0];
                    new_point.y = location[1];
                    new_point.z = location[2];
                    display_rep.points.push_back(new_point);
                    // Color it
                    std_msgs::ColorRGBA new_color;
                    new_color.a = alpha;
                    new_color.b = 0.0;
                    new_color.g = 0.0;
                    new_color.r = 1.0;
                    display_rep.colors.push_back(new_color);
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker SignedDistanceField::ExportForDebug(float alpha) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "sdf_display";
    display_rep.id = 1;
    display_rep.type = visualization_msgs::Marker::CUBE_LIST;
    display_rep.action = visualization_msgs::Marker::ADD;
    display_rep.lifetime = ros::Duration(0.0);
    display_rep.frame_locked = false;
    const Eigen::Affine3d base_transform = Eigen::Affine3d::Identity();
    display_rep.pose = EigenHelpersConversions::EigenAffine3dToGeometryPose(base_transform);
    display_rep.scale.x = GetResolution();
    display_rep.scale.y = GetResolution();
    display_rep.scale.z = GetResolution();
    // Add all the cells of the SDF to the message
    for (int64_t x_index = 0; x_index < distance_field_.GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < distance_field_.GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < distance_field_.GetNumZCells(); z_index++)
            {
                // Convert SDF indices into a real-world location
                std::vector<double> location = distance_field_.GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                display_rep.points.push_back(new_point);
                // Color it
                std_msgs::ColorRGBA new_color;
                new_color.a = alpha;
                new_color.b = 0.0;
                new_color.g = 1.0;
                new_color.r = 1.0;
                display_rep.colors.push_back(new_color);
            }
        }
    }
    return display_rep;
}

void SignedDistanceField::FollowGradientsToLocalMaximaUnsafe(VoxelGrid::VoxelGrid<Eigen::Vector3d>& watershed_map, const int64_t x_index, const int64_t y_index, const int64_t z_index) const
{
    // First, check if we've already found the local maxima for the current cell
    const Eigen::Vector3d& stored = watershed_map.GetImmutable(x_index, y_index, z_index).first;
    if (stored.x() != -INFINITY && stored.y() != -INFINITY && stored.z() != -INFINITY)
    {
        // We've already found it for this cell, so we can skip it
        return;
    }
    // Second, check if it's inside an obstacle
    float stored_distance = Get(x_index, y_index, z_index);
    if (stored_distance <= 0.0)
    {
        // It's inside an object, so we can skip it
        return;
    }
    else
    {
        // Find the local maxima
        std::vector<double> raw_gradient = GetGradient(x_index, y_index, z_index, true);
        Eigen::Vector3d current_gradient(raw_gradient[0], raw_gradient[1], raw_gradient[2]);
        if (GradientIsEffectiveFlat(current_gradient))
        {
            std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
            Eigen::Vector3d local_maxima(location[0], location[1], location[2]);
            watershed_map.SetValue(x_index, y_index, z_index, local_maxima);
        }
        else
        {
            // Follow the gradient, one cell at a time, until we reach a local maxima
            std::unordered_map<VoxelGrid::GRID_INDEX, int8_t> path;
            VoxelGrid::GRID_INDEX current_index(x_index, y_index, z_index);
            path[current_index] = 1;
            Eigen::Vector3d local_maxima(-INFINITY, -INFINITY, -INFINITY);
            while (true)
            {
                if (path.size() == 10000)
                {
                    std::cerr << "Warning, gradient path is long (i.e >= 10000 steps)" << std::endl;
                }
                current_index = GetNextFromGradient(current_index, current_gradient);
                if (path[current_index] != 0)
                {
                    //std::cerr << "LMAX found by cycle detect" << std::endl;
                    // If we've already been here, then we are done
                    std::vector<double> location = GridIndexToLocation(current_index);
                    local_maxima = Eigen::Vector3d(location[0], location[1], location[2]);
                    break;
                }
                // Check if we've been pushed past the edge
                if (current_index.x < 0 || current_index.y < 0 || current_index.z < 0 || current_index.x >= watershed_map.GetNumXCells() || current_index.y >= watershed_map.GetNumYCells() || current_index.z >= watershed_map.GetNumZCells())
                {
                    // We have the "off the grid" local maxima
                    local_maxima = Eigen::Vector3d(INFINITY, INFINITY, INFINITY);
                    break;
                }
                path[current_index] = 1;
                // Check if the new index has already been checked
                const Eigen::Vector3d& new_stored = watershed_map.GetImmutable(current_index).first;
                if (new_stored.x() != -INFINITY && new_stored.y() != -INFINITY && new_stored.z() != -INFINITY)
                {
                    // We have the local maxima
                    local_maxima = new_stored;
                    break;
                }
                else
                {
                    raw_gradient = GetGradient(current_index, true);
                    current_gradient = Eigen::Vector3d(raw_gradient[0], raw_gradient[1], raw_gradient[2]);
                    if (GradientIsEffectiveFlat(current_gradient))
                    {
                        //std::cerr << "LMAX found by flat detect" << std::endl;
                        // We have the local maxima
                        std::vector<double> location = GridIndexToLocation(current_index);
                        local_maxima = Eigen::Vector3d(location[0], location[1], location[2]);
                        break;
                    }
                }
            }
            // Now, go back and mark the entire explored path with the local maxima
            std::unordered_map<VoxelGrid::GRID_INDEX, int8_t>::const_iterator path_itr;
            for (path_itr = path.begin(); path_itr != path.end(); ++path_itr)
            {
                const VoxelGrid::GRID_INDEX& index = path_itr->first;
                watershed_map.SetValue(index, local_maxima);
            }
        }
    }
}

VoxelGrid::VoxelGrid<Eigen::Vector3d> SignedDistanceField::ComputeLocalMaximaMap() const
{
    VoxelGrid::VoxelGrid<Eigen::Vector3d> watershed_map(GetOriginTransform(), GetResolution(), GetXSize(), GetYSize(), GetZSize(), Eigen::Vector3d(-INFINITY, -INFINITY, -INFINITY));
    for (int64_t x_idx = 0; x_idx < watershed_map.GetNumXCells(); x_idx++)
    {
        for (int64_t y_idx = 0; y_idx < watershed_map.GetNumYCells(); y_idx++)
        {
            for (int64_t z_idx = 0; z_idx < watershed_map.GetNumZCells(); z_idx++)
            {
                // We use an "unsafe" function here because we know all the indices we provide it will be safe
                FollowGradientsToLocalMaximaUnsafe(watershed_map, x_idx, y_idx, z_idx);
            }
        }
    }
    return watershed_map;
}
