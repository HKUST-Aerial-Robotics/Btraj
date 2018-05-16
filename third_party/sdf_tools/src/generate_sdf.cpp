#include <ros/ros.h>
#include <arc_utilities/pretty_print.hpp>
#include <arc_utilities/voxel_grid.hpp>
#include "sdf_tools/CollisionMap.h"
#include "sdf_tools/RequestSDF.h"
#include "sdf_tools/SDF.h"
#include "sdf_tools/collision_map.hpp"
#include "sdf_tools/sdf.hpp"

//ssssss

sdf_tools::SignedDistanceField sdf;
// service callback function
bool requestCallback(sdf_tools::RequestSDF::Request &req,
                     sdf_tools::RequestSDF::Response &res)
{
  std::vector<double> location_gradient_query= sdf.GetGradient(req.x, req.y, req.z, true);
  res.gx= location_gradient_query[0];
  res.gy= location_gradient_query[1];
  res.gz= location_gradient_query[2];

  std::pair<float, bool> location_sdf_query= sdf.GetSafe(req.x, req.y, req.z);
  res.distance= location_sdf_query.first;

  return true;
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "generate_sdf");
  ros::NodeHandle nh;

  ros::Publisher visualization_pub= nh.advertise<visualization_msgs::Marker>(
      "sdf_tools_tutorial_visualization", 1, true);

  ros::Publisher sdf_pub= nh.advertise<sdf_tools::SDF>("sdf_pub", 1, true);

  ros::ServiceServer service= nh.advertiseService("request_sdf", requestCallback);

  // wait for advertise
  ros::Duration(1.0).sleep();

  // In preparation, we want to set a couple common paramters
  double resolution= 0.2;
  double x_size= 20.0;
  double y_size= 20.0;
  double z_size= 20.0;
  Eigen::Translation3d origin_translation(-10.0, -10.0, -10.0);
  Eigen::Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
  Eigen::Isometry3d origin_transform= origin_translation * origin_rotation;
  std::string frame= "world";

  // build some obstacle first,start by initialize a collision map
  sdf_tools::COLLISION_CELL oob_cell;
  oob_cell.occupancy= 0.0;
  oob_cell.component= 0;
  sdf_tools::CollisionMapGrid collision_map(origin_transform, frame, resolution, x_size,
                                            y_size, z_size, oob_cell);

  // add collsion cells somewhere
  sdf_tools::COLLISION_CELL obstacle_cell(1.0);

  for(float z= 0; z < 5; z+= resolution)
  {
    // collision_map.Set(-1.5, -0.5, z, obstacle_cell);
    // collision_map.Set(0, 0.5, z, obstacle_cell);
    // collision_map.Set(1.2, 0.3, z, obstacle_cell);
    // collision_map.Set(-0.7, -1.3, z, obstacle_cell);
    // collision_map.Set(0.5, -0.5, z, obstacle_cell);
    collision_map.Set(-1.5, 0.6, z, obstacle_cell);
    collision_map.Set(-0.75, 0.4, z, obstacle_cell);
    collision_map.Set(0, 0.2, z, obstacle_cell);
    collision_map.Set(0.75, 0.4, z, obstacle_cell);
    collision_map.Set(1.5, 0.6, z, obstacle_cell);

    collision_map.Set(-1.5, -0.9, z, obstacle_cell);
    collision_map.Set(-0.75, -1.1, z, obstacle_cell);
    collision_map.Set(0, -1.3, z, obstacle_cell);
    collision_map.Set(0.75, -1.1, z, obstacle_cell);
    collision_map.Set(1.5, -0.9, z, obstacle_cell);
  }

  // visualize collision map in rviz
  std_msgs::ColorRGBA collision_color;
  collision_color.r= 0.0;
  collision_color.g= 0.0;
  collision_color.b= 1.0;
  collision_color.a= 0.8;

  std_msgs::ColorRGBA free_color;
  free_color.r= 0.0;
  free_color.g= 1.0;
  free_color.b= 0.0;
  free_color.a= 0.0;

  std_msgs::ColorRGBA unknown_color;
  unknown_color.r= 1.0;
  unknown_color.g= 1.0;
  unknown_color.b= 0.0;
  unknown_color.a= 0.0;

  visualization_msgs::Marker collision_map_marker=
      collision_map.ExportForDisplay(collision_color, free_color, unknown_color);
  collision_map_marker.ns= "collision_map";
  collision_map_marker.id= 1;

  visualization_pub.publish(collision_map_marker);

  // build a sdf
  float oob_value= INFINITY;
  std::pair<sdf_tools::SignedDistanceField, std::pair<double, double>> sdf_with_extrema=
      collision_map.ExtractSignedDistanceField(oob_value);

  sdf= sdf_with_extrema.first;

  //   // get value from sdf
  //   std::pair<float, bool> location_sdf_query= sdf.GetSafe(-1.5, -0.5, 1);
  //   std::cout << "Location query result - stored distance " << location_sdf_query.first
  //             << " was it in the grid? - " << location_sdf_query.second << std::endl;

  // publish the sdf
  sdf_pub.publish(sdf.GetMessageRepresentation());
  std::cout << "...done" << std::endl;

  ros::spin();

  return 0;
}