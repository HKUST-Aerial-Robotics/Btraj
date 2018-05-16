#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <zlib.h>
#include <ros/ros.h>
#include <list>
#include <unordered_map>
#include <sdf_tools/tagged_object_collision_map.hpp>
#include <arc_utilities/zlib_helpers.hpp>
#include <arc_utilities/eigen_helpers.hpp>
#include <arc_utilities/eigen_helpers_conversions.hpp>
#include <arc_utilities/pretty_print.hpp>
#include <sdf_tools/TaggedObjectCollisionMap.h>

using namespace sdf_tools;

bool TaggedObjectCollisionMapGrid::SaveToFile(const std::string &filepath) const
{
    // Convert to message representation
    sdf_tools::TaggedObjectCollisionMap message_rep = GetMessageRepresentation();
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

bool TaggedObjectCollisionMapGrid::LoadFromFile(const std::string& filepath)
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
        sdf_tools::TaggedObjectCollisionMap new_message;
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

std::vector<uint8_t> TaggedObjectCollisionMapGrid::PackBinaryRepresentation(const std::vector<TAGGED_OBJECT_COLLISION_CELL>& raw) const
{
    std::vector<uint8_t> packed(raw.size() * sizeof(TAGGED_OBJECT_COLLISION_CELL));
    for (size_t field_idx = 0, binary_index = 0; field_idx < raw.size(); field_idx++, binary_index+=sizeof(TAGGED_OBJECT_COLLISION_CELL))
    {
        const TAGGED_OBJECT_COLLISION_CELL& raw_cell = raw[field_idx];
        std::vector<uint8_t> packed_cell = TaggedObjectCollisionCellToBinary(raw_cell);
        memcpy(&packed[binary_index], &packed_cell.front(), sizeof(TAGGED_OBJECT_COLLISION_CELL));
    }
    return packed;
}

std::vector<TAGGED_OBJECT_COLLISION_CELL> TaggedObjectCollisionMapGrid::UnpackBinaryRepresentation(const std::vector<uint8_t>& packed) const
{
    if ((packed.size() % sizeof(TAGGED_OBJECT_COLLISION_CELL)) != 0)
    {
        std::cerr << "Invalid binary representation - length is not a multiple of " << sizeof(TAGGED_OBJECT_COLLISION_CELL) << std::endl;
        return std::vector<TAGGED_OBJECT_COLLISION_CELL>();
    }
    uint64_t data_size = packed.size() / sizeof(TAGGED_OBJECT_COLLISION_CELL);
    std::vector<TAGGED_OBJECT_COLLISION_CELL> unpacked(data_size);
    for (size_t field_idx = 0, binary_index = 0; field_idx < unpacked.size(); field_idx++, binary_index+=sizeof(TAGGED_OBJECT_COLLISION_CELL))
    {
        std::vector<uint8_t> binary_block(sizeof(TAGGED_OBJECT_COLLISION_CELL));
        memcpy(&binary_block.front(), &packed[binary_index], sizeof(TAGGED_OBJECT_COLLISION_CELL));
        unpacked[field_idx] = TaggedObjectCollisionCellFromBinary(binary_block);
    }
    return unpacked;
}

sdf_tools::TaggedObjectCollisionMap TaggedObjectCollisionMapGrid::GetMessageRepresentation() const
{
    sdf_tools::TaggedObjectCollisionMap message_rep;
    // Populate message
    message_rep.header.frame_id = frame_;
    Eigen::Affine3d origin_transform = GetOriginTransform();
    message_rep.origin_transform.translation.x = origin_transform.translation().x();
    message_rep.origin_transform.translation.y = origin_transform.translation().y();
    message_rep.origin_transform.translation.z = origin_transform.translation().z();
    Eigen::Quaterniond origin_transform_rotation(origin_transform.rotation());
    message_rep.origin_transform.rotation.x = origin_transform_rotation.x();
    message_rep.origin_transform.rotation.y = origin_transform_rotation.y();
    message_rep.origin_transform.rotation.z = origin_transform_rotation.z();
    message_rep.origin_transform.rotation.w = origin_transform_rotation.w();
    message_rep.dimensions.x = GetXSize();
    message_rep.dimensions.y = GetYSize();
    message_rep.dimensions.z = GetZSize();
    message_rep.cell_size = GetResolution();
    message_rep.OOB_value = TaggedObjectCollisionCellToBinary(GetOOBValue());
    message_rep.number_of_components = number_of_components_;
    message_rep.components_valid = components_valid_;
    message_rep.convex_segments_valid = convex_segments_valid_;
    message_rep.initialized = initialized_;
    const std::vector<TAGGED_OBJECT_COLLISION_CELL>& raw_data = collision_field_.GetRawData();
    std::vector<uint8_t> binary_data = PackBinaryRepresentation(raw_data);
    message_rep.data = ZlibHelpers::CompressBytes(binary_data);
    return message_rep;
}

bool TaggedObjectCollisionMapGrid::LoadFromMessageRepresentation(const sdf_tools::TaggedObjectCollisionMap& message)
{
    // Make a new voxel grid inside
    Eigen::Translation3d origin_translation(message.origin_transform.translation.x, message.origin_transform.translation.y, message.origin_transform.translation.z);
    Eigen::Quaterniond origin_rotation(message.origin_transform.rotation.w, message.origin_transform.rotation.x, message.origin_transform.rotation.y, message.origin_transform.rotation.z);
    Eigen::Affine3d origin_transform = origin_translation * origin_rotation;
    TAGGED_OBJECT_COLLISION_CELL OOB_value = TaggedObjectCollisionCellFromBinary(message.OOB_value);
    VoxelGrid::VoxelGrid<TAGGED_OBJECT_COLLISION_CELL> new_field(origin_transform, message.cell_size, message.dimensions.x, message.dimensions.y, message.dimensions.z, OOB_value);
    // Unpack the binary data
    std::vector<uint8_t> binary_representation = ZlibHelpers::DecompressBytes(message.data);
    std::vector<TAGGED_OBJECT_COLLISION_CELL> unpacked = UnpackBinaryRepresentation(binary_representation);
    if (unpacked.empty())
    {
        std::cerr << "Unpack returned an empty TaggedObjectCollisionMapGrid" << std::endl;
        return false;
    }
    bool success = new_field.SetRawData(unpacked);
    if (!success)
    {
        std::cerr << "Unable to set internal representation of the TaggedObjectCollisionMapGrid" << std::endl;
        return false;
    }
    // Set it
    collision_field_ = new_field;
    frame_ = message.header.frame_id;
    number_of_components_ = message.number_of_components;
    components_valid_ = message.components_valid;
    convex_segments_valid_ = message.convex_segments_valid;
    initialized_ = message.initialized;
    return true;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportForDisplay(const float alpha, const std::vector<uint32_t>& objects_to_draw) const
{
    std::map<uint32_t, uint32_t> objects_to_draw_map;
    for (size_t idx = 0; idx < objects_to_draw.size(); idx++)
    {
        objects_to_draw_map[objects_to_draw[idx]] = 1u;
    }
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                const auto draw_found_itr = objects_to_draw_map.find(current_cell.object_id);
                if (draw_found_itr != objects_to_draw_map.end() || objects_to_draw_map.size() == 0)
                {
                    const std_msgs::ColorRGBA object_color = GenerateComponentColor(current_cell.object_id, alpha);
                    if (object_color.a > 0.0)
                    {
                        display_rep.points.push_back(new_point);
                        display_rep.colors.push_back(object_color);
                    }
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportForDisplay(const std::map<uint32_t, std_msgs::ColorRGBA>& object_color_map) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                // Check if we've been given a color to work with
                auto found_itr = object_color_map.find(current_cell.object_id);
                std_msgs::ColorRGBA object_color;
                if (found_itr != object_color_map.end())
                {
                    object_color = found_itr->second;
                }
                else
                {
                    object_color = GenerateComponentColor(current_cell.object_id);
                }
                if (object_color.a > 0.0)
                {
                    display_rep.points.push_back(new_point);
                    display_rep.colors.push_back(object_color);
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportContourOnlyForDisplay(const float alpha, const std::vector<uint32_t>& objects_to_draw) const
{
    std::map<uint32_t, uint32_t> objects_to_draw_map;
    for (size_t idx = 0; idx < objects_to_draw.size(); idx++)
    {
        objects_to_draw_map[objects_to_draw[idx]] = 1u;
    }
    // Make SDF
    const std::map<uint32_t, sdf_tools::SignedDistanceField> per_object_sdfs = MakeObjectSDFs();
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                // Get the SDF for the current object
                auto sdf_found_itr = per_object_sdfs.find(current_cell.object_id);
                if (sdf_found_itr != per_object_sdfs.end())
                {
                    const sdf_tools::SignedDistanceField& object_sdf = sdf_found_itr->second;
                    const float distance = object_sdf.Get(new_point.x, new_point.y, new_point.z);
                    // Check if we're on the surface of the object
                    if (distance < 0.0 && distance > -GetResolution())
                    {
                        const auto draw_found_itr = objects_to_draw_map.find(current_cell.object_id);
                        if (draw_found_itr != objects_to_draw_map.end() || objects_to_draw_map.size() == 0)
                        {
                            const std_msgs::ColorRGBA object_color = GenerateComponentColor(current_cell.object_id, alpha);
                            if (object_color.a > 0.0)
                            {
                                display_rep.points.push_back(new_point);
                                display_rep.colors.push_back(object_color);
                            }
                        }
                    }
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportContourOnlyForDisplay(const std::map<uint32_t, std_msgs::ColorRGBA>& object_color_map) const
{
    // Make SDF
    const std::map<uint32_t, sdf_tools::SignedDistanceField> per_object_sdfs = MakeObjectSDFs();
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                // Get the SDF for the current object
                auto sdf_found_itr = per_object_sdfs.find(current_cell.object_id);
                if (sdf_found_itr != per_object_sdfs.end())
                {
                    const sdf_tools::SignedDistanceField& object_sdf = sdf_found_itr->second;
                    const float distance = object_sdf.Get(new_point.x, new_point.y, new_point.z);
                    // Check if we're on the surface of the object
                    if (distance < 0.0 && distance > -GetResolution())
                    {
                        // Check if we've been given a color to work with
                        auto found_itr = object_color_map.find(current_cell.object_id);
                        std_msgs::ColorRGBA object_color;
                        if (found_itr != object_color_map.end())
                        {
                            object_color = found_itr->second;
                        }
                        else
                        {
                            object_color = GenerateComponentColor(current_cell.object_id);
                        }
                        if (object_color.a > 0.0)
                        {
                            display_rep.points.push_back(new_point);
                            display_rep.colors.push_back(object_color);
                        }
                    }
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportForDisplayOccupancyOnly(const std_msgs::ColorRGBA& collision_color, const std_msgs::ColorRGBA& free_color, const std_msgs::ColorRGBA& unknown_color) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_occupancy_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                if (GetImmutable(x_index, y_index, z_index).first.occupancy > 0.5)
                {
                    if (collision_color.a > 0.0)
                    {
                        display_rep.points.push_back(new_point);
                        display_rep.colors.push_back(collision_color);
                    }
                }
                else if (GetImmutable(x_index, y_index, z_index).first.occupancy < 0.5)
                {
                    if (free_color.a > 0.0)
                    {
                        display_rep.points.push_back(new_point);
                        display_rep.colors.push_back(free_color);
                    }
                }
                else
                {
                    if (unknown_color.a > 0.0)
                    {
                        display_rep.points.push_back(new_point);
                        display_rep.colors.push_back(unknown_color);
                    }
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportConnectedComponentsForDisplay(bool color_unknown_components) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_connected_components_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Convert grid indices into a real-world location
                std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                geometry_msgs::Point new_point;
                new_point.x = location[0];
                new_point.y = location[1];
                new_point.z = location[2];
                display_rep.points.push_back(new_point);
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                if (current_cell.occupancy != 0.5)
                {
                    std_msgs::ColorRGBA color = GenerateComponentColor(current_cell.component);
                    display_rep.colors.push_back(color);
                }
                else
                {
                    if (color_unknown_components)
                    {
                        std_msgs::ColorRGBA color = GenerateComponentColor(current_cell.component);
                        display_rep.colors.push_back(color);
                    }
                    else
                    {
                        std_msgs::ColorRGBA color;
                        color.a = 1.0;
                        color.r = 0.5;
                        color.g = 0.5;
                        color.b = 0.5;
                        display_rep.colors.push_back(color);
                    }
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportConvexSegmentForDisplay(const uint32_t object_id, const uint32_t convex_segment) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_" + std::to_string(object_id) + "_convex_segment_" + std::to_string(convex_segment) + "_display";
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
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                const TAGGED_OBJECT_COLLISION_CELL& current_cell = GetImmutable(x_index, y_index, z_index).first;
                if ((current_cell.object_id == object_id) && (current_cell.IsPartOfConvexSegment(convex_segment)))
                {
                    // Convert grid indices into a real-world location
                    std::vector<double> location = GridIndexToLocation(x_index, y_index, z_index);
                    geometry_msgs::Point new_point;
                    new_point.x = location[0];
                    new_point.y = location[1];
                    new_point.z = location[2];
                    display_rep.points.push_back(new_point);
                    // Generate a color
                    const std_msgs::ColorRGBA color = GenerateComponentColor(convex_segment);
                    display_rep.colors.push_back(color);
                }
            }
        }
    }
    return display_rep;
}

visualization_msgs::Marker TaggedObjectCollisionMapGrid::ExportSurfaceForDisplay(const std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>& surface, const std_msgs::ColorRGBA& surface_color) const
{
    // Assemble a visualization_markers::Marker representation of the SDF to display in RViz
    visualization_msgs::Marker display_rep;
    // Populate the header
    display_rep.header.frame_id = frame_;
    // Populate the options
    display_rep.ns = "tagged_object_collision_map_surface";
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
    // Add all the cells of the surface
    std::unordered_map<VoxelGrid::GRID_INDEX, uint8_t>::const_iterator surface_itr;
    for (surface_itr = surface.begin(); surface_itr != surface.end(); ++surface_itr)
    {
        VoxelGrid::GRID_INDEX index = surface_itr->first;
        int8_t validity = surface_itr->second;
        if (validity == 1)
        {
            // Convert grid indices into a real-world location
            std::vector<double> location = GridIndexToLocation(index.x, index.y, index.z);
            geometry_msgs::Point new_point;
            new_point.x = location[0];
            new_point.y = location[1];
            new_point.z = location[2];
            display_rep.points.push_back(new_point);
            display_rep.colors.push_back(surface_color);
        }
    }
    return display_rep;
}

VoxelGrid::VoxelGrid<std::vector<uint32_t>> TaggedObjectCollisionMapGrid::ComputeConvexRegions(const double max_check_radius) const
{
    VoxelGrid::VoxelGrid<std::vector<uint32_t>> convex_region_grid(GetOriginTransform(), GetResolution(), GetXSize(), GetYSize(), GetZSize(), std::vector<uint32_t>());
    uint32_t current_convex_region = 0;
    for (int64_t x_index = 0; x_index < GetNumXCells(); x_index++)
    {
        for (int64_t y_index = 0; y_index < GetNumYCells(); y_index++)
        {
            for (int64_t z_index = 0; z_index < GetNumZCells(); z_index++)
            {
                // Check if cell is empty
                if (GetImmutable(x_index, y_index, z_index).first.occupancy < 0.5)
                {
                    // Check if we've already marked it once
                    const std::vector<uint32_t>& current_cell_regions = convex_region_grid.GetImmutable(x_index, y_index, z_index).first;
                    if (current_cell_regions.empty())
                    {
                        current_convex_region++;
                        std::cout << "Marking convex region " << current_convex_region << std::endl;
                        GrowConvexRegion(VoxelGrid::GRID_INDEX(x_index, y_index, z_index), convex_region_grid, max_check_radius, current_convex_region);
                    }
                }
            }
        }
    }
    std::cout << "Marked " << current_convex_region << " convex regions" << std::endl;
    return convex_region_grid;
}

std::vector<VoxelGrid::GRID_INDEX> TaggedObjectCollisionMapGrid::CheckIfConvex(const VoxelGrid::GRID_INDEX& candidate_index, std::unordered_map<VoxelGrid::GRID_INDEX, int8_t>& explored_indices, const VoxelGrid::VoxelGrid<std::vector<uint32_t>>& region_grid, const uint32_t current_convex_region) const
{
    std::vector<VoxelGrid::GRID_INDEX> convex_indices;
    for (auto indices_itr = explored_indices.begin(); indices_itr != explored_indices.end(); ++indices_itr)
    {
        const VoxelGrid::GRID_INDEX& other_index = indices_itr->first;
        const int8_t& other_status = indices_itr->second;
        // We only care about indices that are already part of the convex set
        if (other_status == 1)
        {
            // Walk from first index to second index. If any intervening cells are filled, return false
            const Eigen::Vector3d start_location = EigenHelpers::StdVectorDoubleToEigenVector3d(GridIndexToLocation(other_index.x, other_index.y, other_index.z));
            const Eigen::Vector3d end_location = EigenHelpers::StdVectorDoubleToEigenVector3d(GridIndexToLocation(candidate_index.x, candidate_index.y, candidate_index.z));
            double distance = (end_location - start_location).norm();
            uint32_t num_steps = (uint32_t)ceil(distance / (GetResolution() * 0.5));
            for (uint32_t step_num = 0; step_num <= num_steps; step_num++)
            {
                const double ratio = (double)step_num / (double)num_steps;
                const Eigen::Vector3d interpolated_location = EigenHelpers::Interpolate(start_location, end_location, ratio);
                std::vector<int64_t> raw_interpolated_index = region_grid.LocationToGridIndex3d(interpolated_location);
                assert(raw_interpolated_index.size() == 3);
                VoxelGrid::GRID_INDEX interpolated_index(raw_interpolated_index[0], raw_interpolated_index[1], raw_interpolated_index[2]);
                // Grab the cell at that location
                const TAGGED_OBJECT_COLLISION_CELL& intermediate_cell = GetImmutable(interpolated_index).first;
                // Check for collision
                if (intermediate_cell.occupancy >= 0.5)
                {
                    return convex_indices;
                }
                // Check if we've already explored it
                if (explored_indices[interpolated_index] == 1)
                {
                    // Great
                    ;
                }
                else if (explored_indices[interpolated_index] == -1)
                {
                    // We've already skipped it deliberately
                    return convex_indices;
                }
                else
                {
                    if (interpolated_index == candidate_index)
                    {
                        // Great
                        ;
                    }
                    else
                    {
                        // We have no idea, let's see if it is convex with our already-explored indices
                        // Temporarily, we mark ourselves as successful
                        explored_indices[candidate_index] = 1;
                        // Call ourselves with the intermediate location
                        std::vector<VoxelGrid::GRID_INDEX> intermediate_convex_indices = CheckIfConvex(interpolated_index, explored_indices, region_grid, current_convex_region);
                        // Unmark ourselves since we don't really know
                        explored_indices[candidate_index] = 0;
                        // Save the intermediate index for addition
                        convex_indices.insert(convex_indices.end(), intermediate_convex_indices.begin(), intermediate_convex_indices.end());
                        // Check if the intermediate index is convex
                        auto is_convex = std::find(convex_indices.begin(), convex_indices.end(), interpolated_index);
                        if (is_convex == convex_indices.end())
                        {
                            // If not, we're done
                            return convex_indices;
                        }
                    }
                }
            }
        }
    }
    // If all indices were reachable, we are part of the convex set
    convex_indices.push_back(candidate_index);
    explored_indices[candidate_index] = 1;
    return convex_indices;
}

void TaggedObjectCollisionMapGrid::GrowConvexRegion(const VoxelGrid::GRID_INDEX& start_index, VoxelGrid::VoxelGrid<std::vector<uint32_t>>& region_grid, const double max_check_radius, const uint32_t current_convex_region) const
{
    // Mark the region of the start index
    region_grid.GetMutable(start_index).first.push_back(current_convex_region);
    std::cout << "Added " << PrettyPrint::PrettyPrint(start_index) << " to region " << current_convex_region << std::endl;
    const Eigen::Vector3d start_location = EigenHelpers::StdVectorDoubleToEigenVector3d(GridIndexToLocation(start_index.x, start_index.y, start_index.z));
    std::list<VoxelGrid::GRID_INDEX> working_queue;
    std::unordered_map<VoxelGrid::GRID_INDEX, int8_t> queued_hashtable;
    working_queue.push_back(start_index);
    queued_hashtable[start_index] = 1;
    while (working_queue.size() > 0)
    {
        // Get the top of the working queue
        VoxelGrid::GRID_INDEX current_index = working_queue.front();
        // Remove from the queue
        working_queue.pop_front();
        // See if we can add the neighbors
        std::vector<VoxelGrid::GRID_INDEX> potential_neighbors(6);
        potential_neighbors[0] = VoxelGrid::GRID_INDEX(current_index.x - 1, current_index.y, current_index.z);
        potential_neighbors[1] = VoxelGrid::GRID_INDEX(current_index.x + 1, current_index.y, current_index.z);
        potential_neighbors[2] = VoxelGrid::GRID_INDEX(current_index.x, current_index.y - 1, current_index.z);
        potential_neighbors[3] = VoxelGrid::GRID_INDEX(current_index.x, current_index.y + 1, current_index.z);
        potential_neighbors[4] = VoxelGrid::GRID_INDEX(current_index.x, current_index.y, current_index.z - 1);
        potential_neighbors[5] = VoxelGrid::GRID_INDEX(current_index.x, current_index.y, current_index.z + 1);
        for (size_t idx = 0; idx < potential_neighbors.size(); idx++)
        {
            const VoxelGrid::GRID_INDEX& candidate_neighbor = potential_neighbors[idx];
            // Make sure the candidate neighbor is in range
            if ((candidate_neighbor.x >= 0) && (candidate_neighbor.y >= 0) && (candidate_neighbor.z >= 0) && (candidate_neighbor.x < GetNumXCells()) && (candidate_neighbor.y < GetNumYCells()) && (candidate_neighbor.z < GetNumZCells()))
            {
                // Make sure it's within the check radius
                const Eigen::Vector3d current_location = EigenHelpers::StdVectorDoubleToEigenVector3d(GridIndexToLocation(candidate_neighbor.x, candidate_neighbor.y, candidate_neighbor.z));
                double distance = (current_location - start_location).norm();
                if (distance <= max_check_radius)
                {
                    // Make sure the candidate neighbor is empty
                    if (GetImmutable(candidate_neighbor).first.occupancy < 0.5)
                    {
                        // Make sure we haven't already checked
                        if (queued_hashtable[candidate_neighbor] == 0)
                        {
                            // Now, let's check if the current index forms a convex set with the indices marked already
                            std::vector<VoxelGrid::GRID_INDEX> convex_indices = CheckIfConvex(candidate_neighbor, queued_hashtable, region_grid, current_convex_region);
                            // Set this to false. If it really is convex, this will get changed in the next loop
                            queued_hashtable[candidate_neighbor] = -1;
                            // Add the new convex indices
                            for (size_t cdx = 0; cdx < convex_indices.size(); cdx++)
                            {
                                const VoxelGrid::GRID_INDEX& convex_index = convex_indices[cdx];
                                // Add to the queue
                                working_queue.push_back(convex_index);
                                queued_hashtable[convex_index] = 1;
                                // Mark it
                                region_grid.GetMutable(convex_index).first.push_back(current_convex_region);
                                std::cout << "Added " << PrettyPrint::PrettyPrint(convex_index) << " to region " << current_convex_region << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
}

