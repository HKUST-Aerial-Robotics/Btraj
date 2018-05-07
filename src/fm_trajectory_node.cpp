#include <iostream>
#include <fstream>
#include <math.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <pcl_conversions/pcl_conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>

#include <ros/ros.h>
#include <ros/console.h>
#include <sensor_msgs/PointCloud2.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>

#include <tf/tf.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>

#include "fm_planer/trajectory_generator.h"
#include "fm_planer/bezier_base.h"
#include "fm_planer/dataType.h"
#include "fm_planer/utils.h"
#include "fm_planer/backward.hpp"

#include "quadrotor_msgs/PositionCommand.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"

using namespace std;
using namespace Eigen;

namespace backward {
backward::SignalHandling sh;
}

// input msgs
nav_msgs::Odometry _odom;
bool _has_odom     = false;
bool _isTargetRcv  = false;

const size_t _odom_queue_size = 200;
deque<nav_msgs::Odometry> _odom_queue;

ros::Subscriber _map_sub, _cmd_sub, _pts_sub, _odom_sub;
ros::Publisher _path_vis_pub, _map_inflation_vis_pub, _corridor_vis_pub, _traj_vis_pub, _traj_pub, _checkTraj_vis_pub, _checkTraj_vis_pcd_pub;
sdf_tools::SignedDistanceField sdf;

double _vis_traj_width = 0.15;
double resolution = 0.2;
double _cloud_margin, _cube_margin;

bool has_path   = true;  // by defalut there is a path
bool has_traj  = false;
bool rcv_sdf    = false;
bool rcv_target = false;

Vector3d startPt, startVel, startAcc, endPt;

Vector3d mapOrigin(-25.0, -25.0, 0.0);
Vector3d local_origin;
double _x_size = 50.0, _y_size = 50.0, _z_size = 4.0;

//unsigned int size_x, size_y, size_z;
unsigned int size_x = (int)_x_size / resolution;
unsigned int size_y = (int)_y_size / resolution;
unsigned int size_z = (int)_z_size / resolution; // fix this when doing experiments
Coord3D dimsize {size_x,size_y,size_z};
FMGrid3D grid_fmm(dimsize);

double pt_max_x = mapOrigin(0) + _x_size;
double pt_min_x = mapOrigin(0);
double pt_max_y = mapOrigin(1) + _y_size;
double pt_min_y = mapOrigin(1); 
double pt_max_z = mapOrigin(2) + _z_size;
double pt_min_z = mapOrigin(2);
    
Translation3d origin_translation( 0.0 + mapOrigin(0), 0.0 + mapOrigin(1), 0.0);
Quaterniond origin_rotation(1.0, 0.0, 0.0, 0.0);
Affine3d origin_transform = origin_translation * origin_rotation;
sdf_tools::COLLISION_CELL oob_cell(0.0);
sdf_tools::CollisionMapGrid collision_map(origin_transform, "world", resolution, _x_size, _y_size, _z_size, oob_cell);
//sdf_tools::CollisionMapGrid collision_map_sdf(origin_transform, "world", resolution, _x_size, _y_size, _z_size, oob_cell);

int max_x = (int)collision_map.GetNumXCells();
int max_y = (int)collision_map.GetNumYCells();
int max_z = (int)collision_map.GetNumZCells();

pcl::PointCloud<pcl::PointXYZ> cloud;

void rcvMapCallBack(sdf_tools::SDF sdf);
void rcvWaypointsCallback(const nav_msgs::Path & wp);
void fastMarching3D();
bool isContains(Cube cube1, Cube cube2);
bool checkHalfWay();
bool checkPointOccupied(vector<double> check_pt);

void visPath(Path3D path);
void visCorridor(vector<Cube> corridor);
void visBezierTrajectory(MatrixXd polyCoeff, VectorXd time);

vector<double> getStateFromBezier(const MatrixXd & polyCoeff, double t_now, int seg_now );
vector<double> getPosFromBezier(const MatrixXd & polyCoeff, double t_now, int seg_now );
void timeAllocation(vector<Cube> & corridor, vector<double> time);
quadrotor_msgs::PolynomialTrajectory getBezierTraj();

TrajectoryGenerator _trajectoryGenerator;
VectorXd _Time;
MatrixXd _PolyCoeff;
MatrixXd _MQM, _FM;
VectorXd _C, _Cv, _Ca, _Cj;
int _segment_num;
int _traj_id = 1;

double _local_rad, _buffer_size, _check_horizon;
double _MAX_Vel, _MAX_Acc;
bool _isLimitVel, _isLimitAcc;
int _step_length, _max_inflate_iter;
int _minimize_order, _traj_order;

quadrotor_msgs::PolynomialTrajectory _traj;

double obj;
ros::Time _start_time = ros::TIME_MAX;

void rcvWaypointsCallback(const nav_msgs::Path & wp)
{     
    if(wp.poses[0].pose.position.z < 0.0)
    return;

    endPt   << wp.poses[0].pose.position.x,
             wp.poses[0].pose.position.y,
             wp.poses[0].pose.position.z;

    rcv_target = true;

    ROS_INFO("[Fast Marching Node] receive the way-points");
    fastMarching3D(); 
}

vector<pcl::PointXYZ> pointInflate( pcl::PointXYZ pt)
{
    int num   = ceil(_cloud_margin / resolution);
    int num_z = num / 2;
    vector<pcl::PointXYZ> infPts;
    pcl::PointXYZ pt_inf;

    //cout<<"inflation num: "<<num<<endl;
    for(int x = -num ; x <= num; x ++ )
        for(int y = -num ; y <= num; y ++ )
            for(int z = -num_z ; z <= num_z; z ++ )
            {
                pt_inf.x = pt.x + x * resolution;
                pt_inf.y = pt.y + y * resolution;
                pt_inf.z = pt.z + z * resolution;

                if( x == num  || y == num  || z ==  num_z 
                ||  x == -num || y == -num || z == -num_z )
                    infPts.push_back( pt_inf );
            }

    return infPts;
}

pcl::PointCloud<pcl::PointXYZ> cloud_inflation;
void rcvPointCloudCallBack(const sensor_msgs::PointCloud2 & pointcloud_map)
{    
    ros::Time time_1 = ros::Time::now();    
    pcl::fromROSMsg(pointcloud_map, cloud);
    
    if((int)cloud.points.size() == 0)
        return;

    sdf_tools::CollisionMapGrid collision_map_global(origin_transform, "world", resolution, _x_size, _y_size, _z_size, oob_cell);

    _local_rad   = 20.0;
    _buffer_size = _MAX_Vel;

    double _x_local_size = _local_rad + _buffer_size;
    double _y_local_size = _local_rad + _buffer_size;
    double _z_local_size = _z_size;

    local_origin << startPt(0) - _x_local_size/2.0, startPt(1) - _y_local_size/2.0, 0.0;

    Translation3d origin_local_translation( local_origin(0), local_origin(1), local_origin(2));
    Quaterniond origin_local_rotation(1.0, 0.0, 0.0, 0.0);

    Affine3d origin_local_transform = origin_local_translation * origin_local_rotation;
    sdf_tools::CollisionMapGrid collision_map_local(origin_local_transform, "world", resolution, _x_local_size, _y_local_size, _z_local_size, oob_cell);

/*    ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time in new a grids map = %f", (time_2 - time_1).toSec());*/

    vector<pcl::PointXYZ> inflatePts;
    pcl::PointCloud<pcl::PointXYZ> cloud_inflation;
    for (int idx = 0; idx < (int)cloud.points.size(); idx++)
    {   
        auto mk = cloud.points[idx];
        pcl::PointXYZ pt(mk.x, mk.y, mk.z);

        if( fabs(pt.x - startPt(0)) > _local_rad / 2.0 || fabs(pt.y - startPt(1)) > _local_rad / 2.0 )
            continue; 
        
        inflatePts = pointInflate(pt);
        //cout<<"inflatePts size: "<<inflatePts.size()<<endl;
        for(int i = 0; i < (int)inflatePts.size(); i++)
        {   
            pcl::PointXYZ inf_pt = inflatePts[i];
            Vector3d addPt(inf_pt.x, inf_pt.y, inf_pt.z);
            sdf_tools::COLLISION_CELL obstacle_cell(1.0); // Occupancy values > 0.5 are obstacles
            collision_map_local.Set3d(addPt, obstacle_cell);
            collision_map_global.Set3d(addPt, obstacle_cell);
            cloud_inflation.push_back(inf_pt);
        }
    }

    cloud_inflation.width = cloud_inflation.points.size();
    cloud_inflation.height = 1;
    cloud_inflation.is_dense = true;
    cloud_inflation.header.frame_id = "world";

    sensor_msgs::PointCloud2 inflateMap;
    pcl::toROSMsg(cloud_inflation, inflateMap);
    _map_inflation_vis_pub.publish(inflateMap);

    collision_map = collision_map_global;

    /*ros::Time time_3 = ros::Time::now();
    ROS_WARN("Time in set up grids map = %f", (time_3 - time_2).toSec());*/

    float oob_value = INFINITY;
    pair<sdf_tools::SignedDistanceField, pair<double, double>> sdf_with_extrema = collision_map_local.ExtractSignedDistanceField(oob_value);
    //pair<sdf_tools::SignedDistanceField, pair<double, double>> sdf_with_extrema = collision_map.ExtractSignedDistanceField(oob_value);
/*    ros::Time time_4 = ros::Time::now();
    ROS_WARN("time in prepare for sdf = %f", (time_4 - time_3).toSec());*/

    sdf = sdf_with_extrema.first;
    rcv_sdf = true;   

    unsigned int idx;
    double max_v = _MAX_Vel;
    vector<unsigned int> obs;            
    Vector3d pt;
    vector<int64_t> pt_idx;
    double occupancy;

    for(unsigned int k = 0; k < size_z; k++)
    {
        for(unsigned int j = 0; j < size_y; j++)
        {
            for(unsigned int i = 0; i < size_x; i++)
            {
                idx = k * size_y * size_x + j * size_x + i;
                pt << i * resolution + mapOrigin(0), 
                      j * resolution + mapOrigin(1), 
                      k * resolution + mapOrigin(2);

                //pt_idx = collision_map_local.LocationToGridIndex(pt(0), pt(1), pt(2)); 
                //if( pt_idx.size() == 3)
                if( fabs(pt(0) - startPt(0)) <= _local_rad / 2.0 && fabs(pt(1) - startPt(1)) <= _local_rad / 2.0)
                    occupancy = sdf.Get( pt(0), pt(1), pt(2));
                else
                    occupancy = max_v;

                occupancy = (occupancy >= max_v) ? max_v : occupancy;
                /*if( k == 0 || k == size_z - 1 || j == 0 || j == size_y - 1 || i == 0 || i == size_x - 1)
                    occupancy = 0.0;*/

                grid_fmm[idx].setOccupancy(occupancy);
                
                if (grid_fmm[idx].isOccupied())
                    obs.push_back(idx);
            }
        }
    }
   
    grid_fmm.setOccupiedCells(std::move(obs));
    grid_fmm.setLeafSize(resolution);

    ros::Time time_5 = ros::Time::now();
    //ROS_WARN("Time consume in set up grid FMM is %f", (time_5 - time_4).toSec() );
    //ROS_WARN("Time consume before planning is %f", (time_5 - time_1).toSec() );
    
    if( has_path == false)
        return; // there is no path exists in the map, unless a new waypoint is sent, no need to try for planning

    if( checkHalfWay() == true )
        fastMarching3D();
    
    ros::Time time_6 = ros::Time::now();
    //ROS_WARN("Time consume in whole pipeline is %f", (time_6 - time_1).toSec() );
    //ROS_BREAK();
}

bool checkHalfWay()
{   
    if(!has_traj) return false;

    //ros::Time time_1 = ros::Time::now();
    vector<double> check_traj_pt;
    vector<double> state;

    visualization_msgs::Marker _traj_vis;

    geometry_msgs::Point pt;
    _traj_vis.header.stamp       = ros::Time::now();
    _traj_vis.header.frame_id    = "world";
    
    pcl::PointCloud<pcl::PointXYZ> checkTraj;
    checkTraj.height = 1;
    checkTraj.is_dense = true;
    
    sensor_msgs::PointCloud2 traj_pcd;

    _traj_vis.ns = "trajectory/trajectory";
    _traj_vis.id = 0;
    _traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    _traj_vis.action = visualization_msgs::Marker::ADD;
    _traj_vis.scale.x = 2.0 * _vis_traj_width;
    _traj_vis.scale.y = 2.0 * _vis_traj_width;
    _traj_vis.scale.z = 2.0 * _vis_traj_width;
    _traj_vis.pose.orientation.x = 0.0;
    _traj_vis.pose.orientation.y = 0.0;
    _traj_vis.pose.orientation.z = 0.0;
    _traj_vis.pose.orientation.w = 1.0;
    _traj_vis.color.r = 0.0;
    _traj_vis.color.g = 0.0;
    _traj_vis.color.b = 1.0;
    _traj_vis.color.a = 2.0;

    double t_s = max(0.0, (_odom.header.stamp - _start_time).toSec());      
    int idx;
    for (idx = 0; idx < _segment_num; ++idx){
      if (t_s > _Time(idx) && idx + 1 < _segment_num)
          t_s -= _Time(idx);
      else break;
    }

    double duration = 0.0;
    double t_ss;
    check_traj_pt.resize(3);
    for(int i = idx; i < _segment_num; i++ )
    {
        t_ss = (i == idx) ? t_s : 0.0;
        for(double t = t_ss; t < _Time(i); t += 0.01){
            double t_d = duration + t - t_ss;
            if( t_d > _check_horizon ) break;
            state = getPosFromBezier( _PolyCoeff, t/_Time(i), i );
            pt.x = check_traj_pt[0] = _Time(i) * state[0]; 
            pt.y = check_traj_pt[1] = _Time(i) * state[1];
            pt.z = check_traj_pt[2] = _Time(i) * state[2];
            
            pcl::PointXYZ pt_point(pt.x, pt.y, pt.z);
            checkTraj.points.push_back(pt_point);
            checkTraj.width = checkTraj.points.size();

            _traj_vis.points.push_back(pt);

            if( checkPointOccupied(check_traj_pt))
            {   
                ROS_ERROR("predicted collision time is %f ahead", t_d);
                _checkTraj_vis_pub.publish(_traj_vis);
                
                pcl::toROSMsg(checkTraj, traj_pcd);
                traj_pcd.header   = _traj_vis.header;
                _checkTraj_vis_pcd_pub.publish(traj_pcd);

                return true;
            }
      }

      duration += _Time(i) - t_ss;
    }

    pcl::toROSMsg(checkTraj, traj_pcd);
    traj_pcd.header  = _traj_vis.header;
    _checkTraj_vis_pcd_pub.publish(traj_pcd);

    _checkTraj_vis_pub.publish(_traj_vis); 

    /*ros::Time time_2 = ros::Time::now();
    ROS_WARN("Time in collision checking is %f", (time_2 - time_1).toSec());*/

    return false;
}

bool checkPointOccupied(vector<double> check_pt)
{
    if(collision_map.Get(check_pt[0], check_pt[1], check_pt[2]).first.occupancy > 0.5)
        return true;
    else
        return false;
}

Vector3i vec2Vec(vector<int64_t> pt_idx)
{
    return Vector3i(pt_idx[0], pt_idx[1], pt_idx[2]);
}

Vector3d vec2Vec(vector<double> pos)
{
    return Vector3d(pos[0], pos[1], pos[2]);
}

pair<Cube, bool> inflate(Cube cube, Cube lstcube)
{   
    Cube cubeMax = cube;

    //cube.printBox();
    // Inflate sequence: left, right, front, back, below, above                                                                                
    MatrixXi vertex_idx(8, 3);
    for (int i = 0; i < 8; i++)
    {   
        double coord_x = max(min(cube.vertex(i, 0), pt_max_x), pt_min_x);
        double coord_y = max(min(cube.vertex(i, 1), pt_max_y), pt_min_y);
        double coord_z = max(min(cube.vertex(i, 2), pt_max_z), pt_min_z);
        Vector3i pt_idx = vec2Vec( collision_map.LocationToGridIndex(coord_x, coord_y, coord_z) );;

        if( collision_map.Get( (int64_t)pt_idx(0), (int64_t)pt_idx(1), (int64_t)pt_idx(2) ).first.occupancy > 0.5 )
            return make_pair(cubeMax, false);
        
        vertex_idx.row(i) = pt_idx;
    }

    // ROS_WARN("assert success");
    int id_x, id_y, id_z;
/*
           P4------------P3 
           /|           /|              ^
          / |          / |              | z
        P1--|---------P2 |              |
         |  P8--------|--p7             |
         | /          | /               /--------> y
         |/           |/               /  
        P5------------P6              / x
*/           
// now is the left side : (p1 -- p4 -- p8 -- p5) face sweep
// ############################################################################################################
    bool loop;

    MatrixXi vertex_idx_lst = vertex_idx;

    /*ROS_WARN("init vertex");
    cout<<"vertex: \n"<<vertex_idx<<endl;*/

    int iter = 0;
    while(iter < _max_inflate_iter)
    {   
        loop  = true; 
        //cout<<"iter: "<<iter<<endl;
        int y_lo = max(0, vertex_idx(0, 1) - _step_length);
        int y_up = min(max_y, vertex_idx(1, 1) + _step_length);

        for(id_y = vertex_idx(0, 1); id_y >= y_lo; id_y-- )
        {   
            if( loop == false) 
                break;
            
            for(id_x = vertex_idx(0, 0); id_x >= vertex_idx(3, 0); id_x-- )
            {    
                if( loop == false) 
                    break;

                for(id_z = vertex_idx(0, 2); id_z >= vertex_idx(4, 2); id_z-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(0, 1) = min(id_y+2, vertex_idx(0, 1));
            vertex_idx(3, 1) = min(id_y+2, vertex_idx(3, 1));
            vertex_idx(7, 1) = min(id_y+2, vertex_idx(7, 1));
            vertex_idx(4, 1) = min(id_y+2, vertex_idx(4, 1));
        }
        else
            vertex_idx(0, 1) = vertex_idx(3, 1) = vertex_idx(7, 1) = vertex_idx(4, 1) = id_y;

        /*ROS_WARN("[corridor inflation] finished the 1st face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/
        
        // now is the right side : (p2 -- p3 -- p7 -- p6) face
        // ############################################################################################################
        loop = true;
        for(id_y = vertex_idx(1, 1); id_y < y_up; id_y++ )
        {   
            if( loop == false) 
                break;
            
            for(id_x = vertex_idx(1, 0); id_x >= vertex_idx(2, 0); id_x-- )
            {
                if( loop == false) 
                    break;

                for(id_z = vertex_idx(1, 2); id_z >= vertex_idx(5, 2); id_z-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(1, 1) = max(id_y-2, vertex_idx(1, 1));
            vertex_idx(2, 1) = max(id_y-2, vertex_idx(2, 1));
            vertex_idx(6, 1) = max(id_y-2, vertex_idx(6, 1));
            vertex_idx(5, 1) = max(id_y-2, vertex_idx(5, 1));
        }
        else
            vertex_idx(1, 1) = vertex_idx(2, 1) = vertex_idx(6, 1) = vertex_idx(5, 1) = id_y;

        /*ROS_WARN("finish the loop of the 2nd face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/

        // now is the front side : (p1 -- p2 -- p6 -- p5) face
        // ############################################################################################################
        int x_lo = max(0, vertex_idx(3, 0) - _step_length);
        int x_up = min(max_x, vertex_idx(0, 0) + _step_length);
        loop = true;
        for(id_x = vertex_idx(0, 0); id_x < x_up; id_x++ )
        {   
            if( loop == false) 
                break;
            
            for(id_y = vertex_idx(0, 1); id_y <= vertex_idx(1, 1); id_y++ )
            {
                if( loop == false) 
                    break;

                for(id_z = vertex_idx(0, 2); id_z >= vertex_idx(4, 2); id_z-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(0, 0) = max(id_x-2, vertex_idx(0, 0)); 
            vertex_idx(1, 0) = max(id_x-2, vertex_idx(1, 0)); 
            vertex_idx(5, 0) = max(id_x-2, vertex_idx(5, 0)); 
            vertex_idx(4, 0) = max(id_x-2, vertex_idx(4, 0)); 
        }
        else
            vertex_idx(0, 0) = vertex_idx(1, 0) = vertex_idx(5, 0) = vertex_idx(4, 0) = id_x;    

        /*ROS_WARN("finish the loop of the 3rd face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/

        // now is the back side : (p4 -- p3 -- p7 -- p8) face
        // ############################################################################################################
        loop = true;
        for(id_x = vertex_idx(3, 0); id_x >= x_lo; id_x-- )
        {   
            if( loop == false) 
                break;
            
            for(id_y = vertex_idx(3, 1); id_y <= vertex_idx(2, 1); id_y++ )
            {
                if( loop == false) 
                    break;

                for(id_z = vertex_idx(3, 2); id_z >= vertex_idx(7, 2); id_z-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(3, 0) = min(id_x+2, vertex_idx(3, 0)); 
            vertex_idx(2, 0) = min(id_x+2, vertex_idx(2, 0)); 
            vertex_idx(6, 0) = min(id_x+2, vertex_idx(6, 0)); 
            vertex_idx(7, 0) = min(id_x+2, vertex_idx(7, 0)); 
        }
        else
            vertex_idx(3, 0) = vertex_idx(2, 0) = vertex_idx(6, 0) = vertex_idx(7, 0) = id_x;

        /*ROS_WARN("finish the loop of the 4th face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/
        // now is the above side : (p1 -- p2 -- p3 -- p4) face
        // ############################################################################################################
        loop = true;
        int z_lo = max(0, vertex_idx(4, 2) - _step_length);
        int z_up = min(max_z, vertex_idx(0, 2) + _step_length);
        for(id_z = vertex_idx(0, 2); id_z < z_up; id_z++ )
        {   
            if( loop == false) 
                break;
            
            for(id_y = vertex_idx(0, 1); id_y <= vertex_idx(1, 1); id_y++ )
            {
                if( loop == false) 
                    break;

                for(id_x = vertex_idx(0, 0); id_x >= vertex_idx(3, 0); id_x-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(0, 2) = max(id_z-2, vertex_idx(0, 2));
            vertex_idx(1, 2) = max(id_z-2, vertex_idx(1, 2));
            vertex_idx(2, 2) = max(id_z-2, vertex_idx(2, 2));
            vertex_idx(3, 2) = max(id_z-2, vertex_idx(3, 2));
        }
        vertex_idx(0, 2) = vertex_idx(1, 2) = vertex_idx(2, 2) = vertex_idx(3, 2) = id_z;

        /*ROS_WARN("finish the loop of the 5th face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/
        // now is the below side : (p5 -- p6 -- p7 -- p8) face
        // ############################################################################################################
        loop = true;
        for(id_z = vertex_idx(4, 2); id_z >= z_lo; id_z-- )
        {   
            if( loop == false) 
                break;
            
            for(id_y = vertex_idx(4, 1); id_y <= vertex_idx(5, 1); id_y++ )
            {
                if( loop == false) 
                    break;

                for(id_x = vertex_idx(4, 0); id_x >= vertex_idx(7, 0); id_x-- )
                {
                    double occupy = collision_map.Get( (int64_t)id_x, (int64_t)id_y, (int64_t)id_z).first.occupancy;    
                    if(occupy > 0.5) // the voxel is occupied
                    {
                        loop = false;
                        break;
                    }
                }
            }
        }

        if(loop == false)
        {
            vertex_idx(4, 2) = min(id_z+2, vertex_idx(4, 2));
            vertex_idx(5, 2) = min(id_z+2, vertex_idx(5, 2));
            vertex_idx(6, 2) = min(id_z+2, vertex_idx(6, 2));
            vertex_idx(7, 2) = min(id_z+2, vertex_idx(7, 2));
        }
        else
            vertex_idx(4, 2) = vertex_idx(5, 2) = vertex_idx(6, 2) = vertex_idx(7, 2) = id_z;

        /*ROS_WARN("finish the loop of the 6th face");
        cout<<"vertex: \n"<<vertex_idx<<endl;*/

        if(vertex_idx_lst == vertex_idx)
        {   
            //ROS_WARN("can not inflate more");
            break;
        }

        vertex_idx_lst = vertex_idx;

        MatrixXd vertex_coord(8, 3);
        for(int i = 0; i < 8; i++)
        {   
            int idx_x = max(min(vertex_idx(i, 0), max_x - 1), 0);
            int idx_y = max(min(vertex_idx(i, 1), max_y - 1), 0);
            int idx_z = max(min(vertex_idx(i, 2), max_z - 1), 0);

            Vector3d pos = vec2Vec( collision_map.GridIndexToLocation( (int64_t)idx_x, (int64_t)idx_y, (int64_t)idx_z) );
            vertex_coord.row(i) = pos;
        }

        cubeMax.setVertex(vertex_coord);
        if(isContains(lstcube, cubeMax))
        {   
            //ROS_WARN("deplicated cube");
            return make_pair(lstcube, false);
        }

        iter ++;
    }

    return make_pair(cubeMax, true);
}

double epsilon = 0.0001;
Cube generateCube( Vector3d pc_) 
{   
/*
           P4------------P3 
           /|           /|              ^
          / |          / |              | z
        P1--|---------P2 |              |
         |  P8--------|--p7             |
         | /          | /               /--------> y
         |/           |/               /  
        P5------------P6              / x
*/       
    Cube cube_;
    
    vector<int64_t> pc_idx    = collision_map.LocationToGridIndex( max(min(pc_(0), pt_max_x), pt_min_x), max(min(pc_(1), pt_max_y), pt_min_y), max(min(pc_(2), pt_max_z), pt_min_z));    
    vector<double>  round_pc_ = collision_map.GridIndexToLocation(pc_idx[0], pc_idx[1], pc_idx[2]);

    cube_.center = Vector3d(round_pc_[0] - resolution/2.0, round_pc_[1] - resolution/2.0, round_pc_[2] - resolution/2.0);

    double x_u = round_pc_[0];//max(min(pc_(0), pt_max_x - epsilon), pt_min_x + epsilon);
    double x_l = round_pc_[0];//max(min(pc_(0), pt_max_x - epsilon), pt_min_x + epsilon);
    
    double y_u = round_pc_[1];//max(min(pc_(1), pt_max_y - epsilon), pt_min_y + epsilon);
    double y_l = round_pc_[1];//max(min(pc_(1), pt_max_y - epsilon), pt_min_y + epsilon);
    
    double z_u = round_pc_[2];//max(min(pc_(2), pt_max_z - epsilon), pt_min_z + epsilon);
    double z_l = round_pc_[2];//max(min(pc_(2), pt_max_z - epsilon), pt_min_z + epsilon);

    cube_.vertex.row(0) = Vector3d(x_u, y_l, z_u);  
    cube_.vertex.row(1) = Vector3d(x_u, y_u, z_u);  
    cube_.vertex.row(2) = Vector3d(x_l, y_u, z_u);  
    cube_.vertex.row(3) = Vector3d(x_l, y_l, z_u);  

    cube_.vertex.row(4) = Vector3d(x_u, y_l, z_l);  
    cube_.vertex.row(5) = Vector3d(x_u, y_u, z_l);  
    cube_.vertex.row(6) = Vector3d(x_l, y_u, z_l);  
    cube_.vertex.row(7) = Vector3d(x_l, y_l, z_l);  

    return cube_;
}

bool isContains(Cube cube1, Cube cube2)
{   
    if( cube1.vertex(0, 0) >= cube2.vertex(0, 0) && cube1.vertex(0, 1) <= cube2.vertex(0, 1) && cube1.vertex(0, 2) >= cube2.vertex(0, 2) &&
        cube1.vertex(6, 0) <= cube2.vertex(6, 0) && cube1.vertex(6, 1) >= cube2.vertex(6, 1) && cube1.vertex(6, 2) <= cube2.vertex(6, 2)  )
        return true;
    else
        return false; 
}

void corridorSimplify(vector<Cube> & cubicList)
{
    vector<Cube> cubicSimplifyList;
    for(int j = (int)cubicList.size() - 1; j > 0; j--)
    {   
        for(int k = j - 1; k > 0; k--)
        {   
            if(cubicList[k].valid == false)
                continue;
            else if(isContains(cubicList[j], cubicList[k]))
                cubicList[k].valid = false;   
        }
    }

    for(auto cube:cubicList)
        if(cube.valid == true)
            cubicSimplifyList.push_back(cube);

    cubicList = cubicSimplifyList;
}

vector<Cube> corridorGeneration(Path3D path, vector<double> time)
{
    vector<Cube> cubeList;
    array<double, 3> state;
    Vector3d pt;

    Cube lstcube;
    lstcube.vertex(0, 0) = -10000;

    for (int i = 0; i < (int)path.size(); i += 1)
    //for (int i = 0; i < 1; i += 1)
    {
        state = path[i];
        pt(0) = state[0];// * resolution + mapOrigin(0);
        pt(1) = state[1];// * resolution + mapOrigin(1);
        pt(2) = state[2];// * resolution;

        if( i > 0 && time[i] == 0.0)
            continue;

        if(isinf(time[i]) || isinf(time[i-1]))
            continue;

        if(i > 0 && time[i] - time[i-1] <= 0.0)
            continue;

        Cube cube = generateCube(pt);
        auto result = inflate(cube, lstcube);
        cube = result.first;
        /*cube.printBox();
        cout<<result.second<<endl;*/

        bool is_skip = false;

        for(int i = 0; i < 3; i ++)
            if( cube.box[i].second - cube.box[i].first < 2 * resolution)
                is_skip = true;

        if( is_skip == true || result.second == false)
            continue;
        
        lstcube = cube;
        cube.t = time[i];
        cubeList.push_back(cube);
    }

    ROS_WARN("Corridor generated, size is %d", (int)cubeList.size() );
    corridorSimplify(cubeList);
    ROS_WARN("Corridor simplified, size is %d", (int)cubeList.size());

    return cubeList;
}

void fastMarching3D()
{   
    if( rcv_target == false || rcv_sdf == false) 
        return;

    Vector3d startIdx3d = (startPt - mapOrigin) / resolution; 
    Vector3d endIdx3d   = (endPt   - mapOrigin) / resolution;

    Coord3D goal_point = {(unsigned int)startIdx3d[0], (unsigned int)startIdx3d[1], (unsigned int)startIdx3d[2]};
    Coord3D init_point = {(unsigned int)endIdx3d[0],   (unsigned int)endIdx3d[1],   (unsigned int)endIdx3d[2]}; 

    unsigned int startIdx;
    vector<unsigned int> startIndices;
    grid_fmm.coord2idx(init_point, startIdx);
    
    startIndices.push_back(startIdx);
    unsigned int goalIdx;
    grid_fmm.coord2idx(goal_point, goalIdx);
    grid_fmm[goalIdx].setOccupancy(0.1); // debug here 
    
    Solver<FMGrid3D>* solver = new FMMStar<FMGrid3D>("FMM*_Dist", TIME); //"FMM*_Dist", DISTANCE
    //Solver<FMGrid3D>* solver = new FMM<FMGrid3D>();

    solver->setEnvironment(&grid_fmm);
    solver->setInitialAndGoalPoints(startIndices, goalIdx);

    ros::Time time_bef_fm = ros::Time::now();
    if(solver->compute(_MAX_Vel) == -1)
    {
        ROS_WARN("[FM Node] No path can be found");
        _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
        _traj_pub.publish(_traj);
        has_traj = false;
        has_path = false;
        return;
    }
    else
        has_path = true;

    Path3D path3D;
    vector<double> path_vels;
    vector<double> time;

    GradientDescent< FMGrid3D > grad3D;
    grid_fmm.coord2idx(goal_point, goalIdx);
    
    double step = 1.0;
    //grad3D.extract_path(grid_fmm, goalIdx, path3D, path_vels, step, time);
    grad3D.apply(grid_fmm, goalIdx, path3D, path_vels, step, time);
    delete solver;

    ros::Time time_aft_fm = ros::Time::now();
    ROS_WARN("[Fast Marching Node] Time in Fast Marching computing is %f", (time_aft_fm - time_bef_fm).toSec() );
    cout << "\tElapsed "<< solver->getName() <<" time: " << solver->getTime() << " ms" << '\n';

    for( int i = 0; i < (int)path3D.size(); i++)
    {
        path3D[i][0] = max(min(path3D[i][0] * resolution + mapOrigin(0), _x_size - resolution), -_x_size + resolution);
        path3D[i][1] = max(min(path3D[i][1] * resolution + mapOrigin(1), _y_size - resolution), -_y_size + resolution);
        path3D[i][2] = max(min(path3D[i][2] * resolution, _z_size - resolution), resolution);
    }
    Vector3d lst3DPt(path3D.back()[0], path3D.back()[1], path3D.back()[2]);
    if((lst3DPt - endPt).norm() > 0.1)
    {
        ROS_WARN("[Fast Marching Node] FMM failed, valid path not exists");
        _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
        _traj_pub.publish(_traj);
        has_traj = false;
        has_path = false;
        return;
    }

    reverse(time.begin(), time.end());
    ros::Time time_bef_corridor = ros::Time::now();
    vector<Cube> corridor = corridorGeneration(path3D, time);

    /*for(auto ptr: corridor)
    {
        ptr.printBox();
        cout<<"time: "<<ptr.t<<endl;
    }*/

    ros::Time time_aft_corridor = ros::Time::now();
    ROS_WARN("Time consume in corridor generation is %f", (time_aft_corridor - time_bef_corridor).toSec());

    MatrixXd pos = MatrixXd::Zero(2,3);
    MatrixXd vel = MatrixXd::Zero(2,3);
    MatrixXd acc = MatrixXd::Zero(2,3);

    pos.row(0) = startPt;
    pos.row(1) = endPt;    
    vel.row(0) = startVel;
    acc.row(0) = startAcc;

    timeAllocation(corridor, time);
    visPath(path3D);
    visCorridor(corridor);

    _segment_num = corridor.size();
    ros::Time time_bef_opt = ros::Time::now();
    /*_PolyCoeff = _trajectoryGenerator.BezierPloyCoeffGeneration(  
                 corridor, _MQM, pos, vel, acc, 2.0, 2.0, _traj_order, _minimize_order, obj, _cube_margin, _isLimitVel, _isLimitAcc );*/
    
    _PolyCoeff = _trajectoryGenerator.BezierPloyCoeffGenerationSOCP(  
                 corridor, _FM, pos, vel, acc, 2.0, 2.0, _traj_order, _minimize_order, obj, _cube_margin, _isLimitVel, _isLimitAcc );
    
    ros::Time time_aft_opt = ros::Time::now();

    ROS_WARN("The objective of the program is %f", obj);
    ROS_WARN("The time consumation of the program is %f", (time_aft_opt - time_bef_opt).toSec());

    if(_PolyCoeff.rows() == 3 && _PolyCoeff.cols() == 3){
          ROS_WARN("Cannot find a feasible and optimal solution, somthing wrong with the mosek solver");
          _traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE;
          _traj_pub.publish(_traj);
          has_traj = false;
          return;
    }
    else
    {
        has_traj = true;
        _traj = getBezierTraj();
        _traj_pub.publish(_traj);
        _traj_id ++;
        ros::Time time_bef_vis = ros::Time::now();
        visBezierTrajectory(_PolyCoeff, _Time);
        ros::Time time_aft_vis = ros::Time::now();
        ROS_WARN("time in visualize the trajectory = %f",(time_aft_vis - time_bef_vis).toSec());
    }
}

void timeAllocation(vector<Cube> & corridor, vector<double> time)
{   
    vector<double> tmp_time;

    for(int i  = 0; i < (int)corridor.size() - 1; i++)
    {   
        tmp_time.push_back(corridor[i+1].t - corridor[i].t);
    }    
    double lst_time  = time.back() - corridor.back().t;

    tmp_time.push_back(lst_time);

/*    cout<<"check time"<<endl;
    for(auto ptr:time)
        cout<<ptr<<endl;
*/
    cout<<"check tmp_time"<<endl;
    for(auto ptr:tmp_time)
        cout<<ptr<<endl;

    _Time.resize((int)corridor.size());

    for(int i = 0; i < (int)corridor.size(); i++)
        _Time(i) = corridor[i].t = tmp_time[i];

    cout<<"allocated time:\n"<<_Time<<endl;
}

void rcvPosCmdCallBack(const quadrotor_msgs::PositionCommand cmd)
{
    startAcc(0)  = cmd.acceleration.x;
    startAcc(1)  = cmd.acceleration.y;
    startAcc(2)  = cmd.acceleration.z;
}

void rcvOdometryCallbck(const nav_msgs::Odometry odom)
{
    if (odom.child_frame_id == "X" || odom.child_frame_id == "O") return ;
    
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    transform.setOrigin( tf::Vector3(_odom.pose.pose.position.x, _odom.pose.pose.position.y, _odom.pose.pose.position.z) );
    transform.setRotation(tf::Quaternion(0, 0, 0, 1.0));

    int i = 0;
    while(i < 10)
    {
        br.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "quadrotor"));
        i++;
    }
    
    _odom = odom;
    _has_odom = true;
    startPt(0)  = _odom.pose.pose.position.x;
    startPt(1)  = _odom.pose.pose.position.y;
    startPt(2)  = _odom.pose.pose.position.z;    

    startVel(0)  = _odom.twist.twist.linear.x;
    startVel(1)  = _odom.twist.twist.linear.y;
    startVel(2)  = _odom.twist.twist.linear.z;    

    _odom_queue.push_back(odom);
    while (_odom_queue.size() > _odom_queue_size) _odom_queue.pop_front();
}

int main(int argc, char** argv)
{
    ros::init(argc, argv, "fast_marching_node");
    ros::NodeHandle nh("~");

    _map_sub   = nh.subscribe( "map", 1,  rcvPointCloudCallBack );
    _cmd_sub   = nh.subscribe( "command", 1,  rcvPosCmdCallBack );
    _odom_sub  = nh.subscribe( "odometry", 50, rcvOdometryCallbck);
    _pts_sub   = nh.subscribe( "waypoints",  1, rcvWaypointsCallback );

    _path_vis_pub          = nh.advertise<visualization_msgs::Marker>("path_vis", 1);
    _map_inflation_vis_pub = nh.advertise<sensor_msgs::PointCloud2>("vis_map_inflate", 1);
    _traj_vis_pub          = nh.advertise<visualization_msgs::Marker>("trajectory_vis", 1);    
    _corridor_vis_pub      = nh.advertise<visualization_msgs::MarkerArray>("corridor_vis", 1);
    _checkTraj_vis_pub     = nh.advertise<visualization_msgs::Marker>("check_trajectory", 10);
    _checkTraj_vis_pcd_pub = nh.advertise<sensor_msgs::PointCloud2>("check_trajectory_pcd", 10);

    _traj_pub = nh.advertise<quadrotor_msgs::PolynomialTrajectory>("trajectory", 10);

    nh.param("optimization/minimize_order", _minimize_order, 3);
    nh.param("optimization/poly_order",     _traj_order,    10);
    nh.param("planning/max_vel",            _MAX_Vel,  1.0);
    nh.param("planning/max_acc",            _MAX_Acc,  1.0);
    nh.param("map/margin",                  _cloud_margin, 0.25);
    nh.param("planning/max_inflate_iter",   _max_inflate_iter, 100);
    nh.param("planning/step_length",        _step_length, 2);
    nh.param("planning/cube_margin",        _cube_margin, 0.2);
    nh.param("planning/check_horizon",      _check_horizon, 10.0);
    nh.param("planning/isLimitVel",         _isLimitVel, false);
    nh.param("planning/isLimitAcc",         _isLimitAcc, false);

    Bernstein _bernstein;
    if(_bernstein.setParam(3, 12, _minimize_order) == -1)
        ROS_ERROR(" The trajectory order is set beyond the library's scope, please re-set "); 

    _MQM = _bernstein.getMQM()[_traj_order];
    _FM  = _bernstein.getFM()[_traj_order];
    _C   = _bernstein.getC()[_traj_order];
    _Cv  = _bernstein.getC_v()[_traj_order];
    _Ca  = _bernstein.getC_a()[_traj_order];
    _Cj  = _bernstein.getC_j()[_traj_order];

    ros::Rate rate(100);
    bool status = ros::ok();
    while(status) 
    {
        ros::spinOnce();           
        status = ros::ok();
        rate.sleep();
    }

    return 0;
}

quadrotor_msgs::PolynomialTrajectory getBezierTraj()
{
    quadrotor_msgs::PolynomialTrajectory traj;
      traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;
      //int _segment_num = _PolyCoeff.rows();
      traj.num_segment = _segment_num;

      int order = _traj_order;
      int poly_num1d = order + 1;
      int polyTotalNum = _segment_num * (order + 1);

      traj.coef_x.resize(polyTotalNum);
      traj.coef_y.resize(polyTotalNum);
      traj.coef_z.resize(polyTotalNum);

      int idx = 0;
      for(int i = 0; i < _segment_num; i++ )
      {    
          for(int j =0; j < poly_num1d; j++)
          { 
              traj.coef_x[idx] = _PolyCoeff(i,                  j);
              traj.coef_y[idx] = _PolyCoeff(i,     poly_num1d + j);
              traj.coef_z[idx] = _PolyCoeff(i, 2 * poly_num1d + j);
              idx++;
          }
      }

      traj.header.frame_id = "/bernstein";
      traj.header.stamp = _odom.header.stamp; //ros::Time(_odom.header.stamp.toSec()); 
      _start_time = traj.header.stamp;

      traj.time.resize(_segment_num);
      traj.order.resize(_segment_num);

      traj.mag_coeff = 1.0;
      for (int idx = 0; idx < _segment_num; ++idx){
          traj.time[idx] = _Time(idx);
          traj.order[idx] = _traj_order;
      }
      
      traj.start_yaw = 0.0;
      traj.final_yaw = 0.0;

      traj.trajectory_id = _traj_id;
      traj.action = quadrotor_msgs::PolynomialTrajectory::ACTION_ADD;

      return traj;
}

vector<double> getPosFromBezier(const MatrixXd & polyCoeff, double t_now, int seg_now )
{
    vector<double > ret(3, 0);
    VectorXd ctrl_now = polyCoeff.row(seg_now);
    int ctrl_num1D = polyCoeff.cols() / 3;

    for(int i = 0; i < 3; i++)
        for(int j = 0; j < ctrl_num1D; j++)
            ret[i] += _C(j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (_traj_order - j) ); 

    return ret;  
}

vector<double> getStateFromBezier(const MatrixXd & polyCoeff, double t_now, int seg_now )
{
    vector<double > ret(12, 0);
    VectorXd ctrl_now = polyCoeff.row(seg_now);
    int ctrl_num1D = polyCoeff.cols() / 3;

    for(int i = 0; i < 3; i++)
    {   
        for(int j = 0; j < ctrl_num1D; j++){
            ret[i] += _C(j) * ctrl_now(i * ctrl_num1D + j) * pow(t_now, j) * pow((1 - t_now), (_traj_order - j) ); 
          
            if(j < ctrl_num1D - 1 )
                ret[i+3] += _Cv(j) * _traj_order 
                      * ( ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j))
                      * pow(t_now, j) * pow((1 - t_now), (_traj_order - j - 1) ); 
          
            if(j < ctrl_num1D - 2 )
                ret[i+6] += _Ca(j) * _traj_order * (_traj_order - 1) 
                      * ( ctrl_now(i * ctrl_num1D + j + 2) - 2 * ctrl_now(i * ctrl_num1D + j + 1) + ctrl_now(i * ctrl_num1D + j))
                      * pow(t_now, j) * pow((1 - t_now), (_traj_order - j - 2) );                         

            if(j < ctrl_num1D - 3 )
                ret[i+9] += _Cj(j) * _traj_order * (_traj_order - 1) * (_traj_order - 2) 
                      * ( ctrl_now(i * ctrl_num1D + j + 3) - 3 * ctrl_now(i * ctrl_num1D + j + 2) + 3 * ctrl_now(i * ctrl_num1D + j + 1) - ctrl_now(i * ctrl_num1D + j))
                      * pow(t_now, j) * pow((1 - t_now), (_traj_order - j - 3) );                         
        }
    }

    return ret;  
}

void visPath(Path3D path)
{
    visualization_msgs::Marker path_vis;
    path_vis.header.stamp       = ros::Time::now();
    path_vis.header.frame_id    = "world";

    path_vis.ns = "trajectory/trajectory";
    path_vis.id = 0;
    path_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    path_vis.action = visualization_msgs::Marker::ADD;
    path_vis.scale.x = _vis_traj_width;
    path_vis.scale.y = _vis_traj_width;
    path_vis.scale.z = _vis_traj_width;
    path_vis.pose.orientation.x = 0.0;
    path_vis.pose.orientation.y = 0.0;
    path_vis.pose.orientation.z = 0.0;
    path_vis.pose.orientation.w = 1.0;
    path_vis.color.r = 1.0;
    path_vis.color.g = 1.0;
    path_vis.color.b = 1.0;
    path_vis.color.a = 1.0;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    path_vis.points.clear();

    array<double, 3> state;
    geometry_msgs::Point pt;

    for (int i = 0; i < (int)path.size(); i += 1, count += 1){
        state = path[i];
        cur(0) = pt.x = state[0];// * resolution + mapOrigin(0);
        cur(1) = pt.y = state[1];// * resolution + mapOrigin(1);
        cur(2) = pt.z = state[2];// * resolution;
        path_vis.points.push_back(pt);

        if (count) traj_len += (pre - cur).norm();
        pre = cur;
    }

    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
    _path_vis_pub.publish(path_vis);
}

visualization_msgs::MarkerArray cube_vis;
void visCorridor(vector<Cube> corridor)
{   
    for(auto & mk: cube_vis.markers) 
        mk.action = visualization_msgs::Marker::DELETE;
    
    _corridor_vis_pub.publish(cube_vis);

    cube_vis.markers.clear();

    visualization_msgs::Marker mk;
    mk.header.frame_id = "world";
    mk.header.stamp = ros::Time::now();
    mk.ns = "corridor";
    mk.type = visualization_msgs::Marker::CUBE;
    mk.action = visualization_msgs::Marker::ADD;

    mk.pose.orientation.x = 0.0;
    mk.pose.orientation.y = 0.0;
    mk.pose.orientation.z = 0.0;
    mk.pose.orientation.w = 1.0;

    mk.color.a = 0.7;
    mk.color.r = 1.0;
    mk.color.g = 1.0;
    mk.color.b = 1.0;

    int idx = 0;
    for(int i = 0; i < int(corridor.size()); i++)
    {   
        mk.id = idx;

        mk.pose.position.x = (corridor[i].vertex(0, 0) + corridor[i].vertex(3, 0) ) / 2.0; 
        mk.pose.position.y = (corridor[i].vertex(0, 1) + corridor[i].vertex(1, 1) ) / 2.0; 
        mk.pose.position.z = -0.1;//(corridor[i].vertex(0, 2) + corridor[i].vertex(4, 2) ) / 2.0; 

        mk.scale.x = (corridor[i].vertex(0, 0) - corridor[i].vertex(3, 0) );
        mk.scale.y = (corridor[i].vertex(1, 1) - corridor[i].vertex(0, 1) );
        mk.scale.z = 0.1;//(corridor[i].vertex(0, 2) - corridor[i].vertex(4, 2) );

        idx ++;
        cube_vis.markers.push_back(mk);
    }

    _corridor_vis_pub.publish(cube_vis);
}

void visBezierTrajectory(MatrixXd polyCoeff, VectorXd time)
{   
    visualization_msgs::Marker traj_vis;

    traj_vis.header.stamp       = ros::Time::now();
    traj_vis.header.frame_id    = "world";

    traj_vis.ns = "trajectory/trajectory";
    traj_vis.id = 0;
    traj_vis.type = visualization_msgs::Marker::SPHERE_LIST;
    
    traj_vis.action = visualization_msgs::Marker::DELETE;
    _checkTraj_vis_pub.publish(traj_vis);

    traj_vis.action = visualization_msgs::Marker::ADD;
    traj_vis.scale.x = _vis_traj_width;
    traj_vis.scale.y = _vis_traj_width;
    traj_vis.scale.z = _vis_traj_width;
    traj_vis.pose.orientation.x = 0.0;
    traj_vis.pose.orientation.y = 0.0;
    traj_vis.pose.orientation.z = 0.0;
    traj_vis.pose.orientation.w = 1.0;
    traj_vis.color.r = 1.0;
    traj_vis.color.g = 0.0;
    traj_vis.color.b = 0.0;
    traj_vis.color.a = 0.6;

    double traj_len = 0.0;
    int count = 0;
    Vector3d cur, pre;
    cur.setZero();
    pre.setZero();
    
    traj_vis.points.clear();

    vector<double> state;
    geometry_msgs::Point pt;

    int segment_num  = polyCoeff.rows();
    for(int i = 0; i < segment_num; i++ ){
        for (double t = 0.0; t < 1.0; t += 0.05 / time(i), count += 1){
            state = getPosFromBezier( polyCoeff, t, i );
            cur(0) = pt.x = time(i) * state[0];
            cur(1) = pt.y = time(i) * state[1];
            cur(2) = pt.z = time(i) * state[2];
            traj_vis.points.push_back(pt);

            if (count) traj_len += (pre - cur).norm();
            pre = cur;
        }
    }

    ROS_INFO("[GENERATOR] The length of the trajectory; %.3lfm.", traj_len);
    _traj_vis_pub.publish(traj_vis);
}