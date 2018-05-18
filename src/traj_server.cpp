#include "ros/ros.h"
#include "ros/console.h"
#include "nav_msgs/Odometry.h"
#include "quadrotor_msgs/PolynomialTrajectory.h"
#include "quadrotor_msgs/PositionCommand.h"
#include "geometry_msgs/PoseStamped.h"
#include "tf/tf.h"
#include "tf/transform_datatypes.h"
#include "bezier_base.h"
#include <eigen3/Eigen/Dense>
#include "visualization_msgs/MarkerArray.h"
#include "visualization_msgs/Marker.h"

const int  _DIM_x = 0;
const int  _DIM_y = 1;
const int  _DIM_z = 2;

using namespace std;

int _poly_order_min, _poly_order_max;

class TrajectoryServer
{
private:

    // Subscribers
    ros::Subscriber _odom_sub;
    ros::Subscriber _traj_sub;

    // publishers
    ros::Publisher _cmd_pub;
    ros::Publisher _vis_cmd_pub;
    ros::Publisher _vis_vel_pub;
    ros::Publisher _vis_acc_pub;
    
    // configuration for trajectory
    int _n_segment = 0;
    int _traj_id = 0;
    uint32_t _traj_flag = 0;
    Eigen::VectorXd _time;
    Eigen::MatrixXd _coef[3];
    vector<int> _order;

    double mag_coeff;
    ros::Time _final_time = ros::TIME_MIN;
    ros::Time _start_time = ros::TIME_MAX;
    double _start_yaw = 0.0, _final_yaw = 0.0;

    // state of the server
    enum ServerState{INIT, TRAJ, HOVER} state = INIT;;
    nav_msgs::Odometry _odom;
    quadrotor_msgs::PositionCommand _cmd;
    geometry_msgs::PoseStamped _vis_cmd;

    visualization_msgs::Marker _vis_vel, _vis_acc;
    Eigen::Matrix3d rot_yaw;

    double vel_t = 0.0;

public:
    
    vector<Eigen::VectorXd> CList;  // Position coefficients vector, used to record all the pre-compute 'n choose k' combinatorial for the bernstein coefficients .
    vector<Eigen::VectorXd> CvList; // Velocity coefficients vector.
    vector<Eigen::VectorXd> CaList; // Acceleration coefficients vector.

    TrajectoryServer(ros::NodeHandle & handle)
    {   
        _odom_sub = 
            handle.subscribe("odometry", 50, &TrajectoryServer::rcvOdometryCallback, this);

        _traj_sub =
            handle.subscribe("trajectory", 2, &TrajectoryServer::rcvTrajectoryCallabck, this);

        _cmd_pub = 
            handle.advertise<quadrotor_msgs::PositionCommand>("position_command", 50);

        _vis_cmd_pub = 
            handle.advertise<geometry_msgs::PoseStamped>("desired_position", 50);

        _vis_vel_pub = 
            handle.advertise<visualization_msgs::Marker>("desired_velocity", 50);
        
        _vis_acc_pub = 
            handle.advertise<visualization_msgs::Marker>("desired_acceleration", 50);

        double pos_gain[3] = {5.7, 5.7, 6.2};
        double vel_gain[3] = {3.4, 3.4, 4.0};
        setGains(pos_gain, vel_gain);
    }

    void setGains(double pos_gain[3], double vel_gain[3])
    {
        _cmd.kx[_DIM_x] = pos_gain[_DIM_x];
        _cmd.kx[_DIM_y] = pos_gain[_DIM_y];
        _cmd.kx[_DIM_z] = pos_gain[_DIM_z];

        _cmd.kv[_DIM_x] = vel_gain[_DIM_x];
        _cmd.kv[_DIM_y] = vel_gain[_DIM_y];
        _cmd.kv[_DIM_z] = vel_gain[_DIM_z];
    }

    void rcvOdometryCallback(const nav_msgs::Odometry & odom)
    {
        if (odom.child_frame_id == "X" || odom.child_frame_id == "O") return ;
        // #1. store the odometry
        _odom = odom;
        _vis_cmd.header = _odom.header;
        _vis_cmd.header.frame_id = "/world";

        if(state == INIT )
        {
            //ROS_WARN("[TRAJ SERVER] Pub initial pos command");
            _cmd.position   = _odom.pose.pose.position;
            
            _cmd.header.stamp = _odom.header.stamp;
            _cmd.header.frame_id = "/world";
            _cmd.trajectory_flag = _traj_flag;

            _cmd.velocity.x = 0.0;
            _cmd.velocity.y = 0.0;
            _cmd.velocity.z = 0.0;
            
            _cmd.acceleration.x = 0.0;
            _cmd.acceleration.y = 0.0;
            _cmd.acceleration.z = 0.0;
            _cmd_pub.publish(_cmd);

            _vis_cmd.pose.position.x = _cmd.position.x;
            _vis_cmd.pose.position.y = _cmd.position.y;
            _vis_cmd.pose.position.z = _cmd.position.z;
            _vis_cmd_pub.publish(_vis_cmd);

            return;
        }
        // #2. try to publish command
        pubPositionCommand();

        // #3. try to calculate the new state
        if (state == TRAJ && ( (odom.header.stamp - _start_time).toSec() / mag_coeff > (_final_time - _start_time).toSec() ) )
        {
            state = HOVER;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
        }
    }

    void rcvTrajectoryCallabck(const quadrotor_msgs::PolynomialTrajectory & traj)
    {
        if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_ADD)
        {   
            ROS_WARN("[SERVER] Loading the trajectory.");
            if ((int)traj.trajectory_id < _traj_id) return ;

            state = TRAJ;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_READY;
            _traj_id = traj.trajectory_id;
            _n_segment = traj.num_segment;
            _final_time = _start_time = traj.header.stamp;
            _time.resize(_n_segment);

            _order.clear();
            for (int idx = 0; idx < _n_segment; ++idx)
            {
                _final_time += ros::Duration(traj.time[idx]);
                _time(idx) = traj.time[idx];
                _order.push_back(traj.order[idx]);
            }

            _start_yaw = traj.start_yaw;
            _final_yaw = traj.final_yaw;
            mag_coeff  = traj.mag_coeff;

            int max_order = *max_element( begin( _order ), end( _order ) ); 
            
            _coef[_DIM_x] = MatrixXd::Zero(max_order + 1, _n_segment);
            _coef[_DIM_y] = MatrixXd::Zero(max_order + 1, _n_segment);
            _coef[_DIM_z] = MatrixXd::Zero(max_order + 1, _n_segment);
            
            int shift = 0;
            for (int idx = 0; idx < _n_segment; ++idx)
            {     
                int order = traj.order[idx];

                for (int j = 0; j < (order + 1); ++j)
                {
                    _coef[_DIM_x](j, idx) = traj.coef_x[shift + j];
                    _coef[_DIM_y](j, idx) = traj.coef_y[shift + j];
                    _coef[_DIM_z](j, idx) = traj.coef_z[shift + j];
                }

                shift += (order + 1);
            }
        }
        else if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_ABORT) 
        {
            ROS_WARN("[SERVER] Aborting the trajectory.");
            state = HOVER;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_COMPLETED;
        }
        else if (traj.action == quadrotor_msgs::PolynomialTrajectory::ACTION_WARN_IMPOSSIBLE)
        {
            state = HOVER;
            _traj_flag = quadrotor_msgs::PositionCommand::TRAJECTORY_STATUS_IMPOSSIBLE;
        }
    }

    void pubPositionCommand()
    {
        if (state == INIT) return;
        if (state == HOVER)
        {
            if (_cmd.header.frame_id != "/world"){
                _cmd.position = _odom.pose.pose.position;
            }

            _cmd.header.stamp = _odom.header.stamp;
            _cmd.header.frame_id = "/world";
            _cmd.trajectory_flag = _traj_flag;

            _cmd.velocity.x = 0.0;
            _cmd.velocity.y = 0.0;
            _cmd.velocity.z = 0.0;
            
            _cmd.acceleration.x = 0.0;
            _cmd.acceleration.y = 0.0;
            _cmd.acceleration.z = 0.0;
        }

        if (state == TRAJ)
        {
            _cmd.header.stamp = _odom.header.stamp;

            _cmd.header.frame_id = "/world";
            _cmd.trajectory_flag = _traj_flag;
            _cmd.trajectory_id = _traj_id;

            double t = max(0.0, (_odom.header.stamp - _start_time).toSec());

            _cmd.yaw_dot = 0.0;
            _cmd.yaw = _start_yaw + (_final_yaw - _start_yaw) * t / ((_final_time - _start_time).toSec() + 1e-9);

            for (int idx = 0; idx < _n_segment; ++idx)
            {
                if (t > _time[idx] && idx + 1 < _n_segment)
                {
                    t -= _time[idx];
                }
                else
                {   
                        t /= _time[idx];

                        _cmd.position.x = 0.0;
                        _cmd.position.y = 0.0;
                        _cmd.position.z = 0.0;
                        _cmd.velocity.x = 0.0;
                        _cmd.velocity.y = 0.0;
                        _cmd.velocity.z = 0.0;
                        _cmd.acceleration.x = 0.0;
                        _cmd.acceleration.y = 0.0;
                        _cmd.acceleration.z = 0.0;

                        int cur_order = _order[idx];
                        int cur_poly_num = cur_order + 1;

                        for(int i = 0; i < cur_poly_num; i ++)
                        {
                            _cmd.position.x += _time[idx] * CList[cur_order](i) * _coef[_DIM_x].col(idx)(i) * pow(t, i) * pow((1 - t), (cur_order - i) ); 
                            _cmd.position.y += _time[idx] * CList[cur_order](i) * _coef[_DIM_y].col(idx)(i) * pow(t, i) * pow((1 - t), (cur_order - i) ); 
                            _cmd.position.z += _time[idx] * CList[cur_order](i) * _coef[_DIM_z].col(idx)(i) * pow(t, i) * pow((1 - t), (cur_order - i) ); 

                            if(i < (cur_poly_num - 1))
                            {
                                _cmd.velocity.x += CvList[cur_order](i) * cur_order * (_coef[_DIM_x].col(idx)(i+1) - _coef[_DIM_x].col(idx)(i)) 
                                                * pow(t, i) * pow((1 - t), (cur_order - 1 - i) );
                                
                                _cmd.velocity.y += CvList[cur_order](i) * cur_order * (_coef[_DIM_y].col(idx)(i+1) - _coef[_DIM_y].col(idx)(i)) 
                                                * pow(t, i) * pow((1 - t), (cur_order - 1 - i) ); 
                                
                                _cmd.velocity.z += CvList[cur_order](i) * cur_order * (_coef[_DIM_z].col(idx)(i+1) - _coef[_DIM_z].col(idx)(i)) 
                                                * pow(t, i) * pow((1 - t), (cur_order - 1 - i) ); 
                            }

                            if(i < (cur_poly_num - 2))
                            {   
                                _cmd.acceleration.x += 1.0 / _time[idx] * CaList[cur_order](i) * cur_order * (cur_order - 1) 
                                                    * (_coef[_DIM_x].col(idx)(i+2) - 2*_coef[_DIM_x].col(idx)(i+1) + _coef[_DIM_x].col(idx)(i)) 
                                                    * pow(t, i) * pow((1 - t), (cur_order - 2 - i) );
                                
                                _cmd.acceleration.y += 1.0 / _time[idx] * CaList[cur_order](i) * cur_order * (cur_order - 1) 
                                                    * (_coef[_DIM_y].col(idx)(i+2) - 2*_coef[_DIM_y].col(idx)(i+1) + _coef[_DIM_y].col(idx)(i)) 
                                                    * pow(t, i) * pow((1 - t), (cur_order - 2 - i) ); 
                                
                                _cmd.acceleration.z += 1.0 / _time[idx] * CaList[cur_order](i) * cur_order * (cur_order - 1) 
                                                    * (_coef[_DIM_z].col(idx)(i+2) - 2*_coef[_DIM_z].col(idx)(i+1) + _coef[_DIM_z].col(idx)(i)) 
                                                    * pow(t, i) * pow((1 - t), (cur_order - 2 - i) ); 
                            }

                        }
                    break;
                } 
            }
        }

        _cmd_pub.publish(_cmd);

        _vis_cmd.header = _cmd.header;
        _vis_cmd.pose.position.x = _cmd.position.x;
        _vis_cmd.pose.position.y = _cmd.position.y;
        _vis_cmd.pose.position.z = _cmd.position.z;
        
        tf::Quaternion q_ = tf::createQuaternionFromYaw(_cmd.yaw);
        geometry_msgs::Quaternion odom_quat;
        tf::quaternionTFToMsg(q_, odom_quat);
        _vis_cmd.pose.orientation = odom_quat;
        _vis_cmd_pub.publish(_vis_cmd);

        _vis_vel.ns = "vel";
        _vis_vel.id = 0;
        _vis_vel.header.frame_id = "/world";
        _vis_vel.type = visualization_msgs::Marker::ARROW;
        _vis_vel.action = visualization_msgs::Marker::ADD;
        _vis_vel.color.a = 1.0;
        _vis_vel.color.r = 0.0;
        _vis_vel.color.g = 1.0;
        _vis_vel.color.b = 0.0;

        _vis_vel.header.stamp = _odom.header.stamp;

        _vis_vel.points.clear();
        geometry_msgs::Point pt;
        pt.x = _cmd.position.x;
        pt.y = _cmd.position.y;
        pt.z = _cmd.position.z;

        _vis_vel.points.push_back(pt);
        
        pt.x = _cmd.position.x + _cmd.velocity.x;
        pt.y = _cmd.position.y + _cmd.velocity.y;
        pt.z = _cmd.position.z + _cmd.velocity.z;

        _vis_vel.points.push_back(pt);

        _vis_vel.scale.x = 0.1;
        _vis_vel.scale.y = 0.2;
        _vis_vel.scale.z = 0.2;

        _vis_vel_pub.publish(_vis_vel);

        _vis_acc.ns = "acc";
        _vis_acc.id = 0;
        _vis_acc.header.frame_id = "/world";
        _vis_acc.type = visualization_msgs::Marker::ARROW;
        _vis_acc.action = visualization_msgs::Marker::ADD;
        _vis_acc.color.a = 1.0;
        _vis_acc.color.r = 1.0;
        _vis_acc.color.g = 1.0;
        _vis_acc.color.b = 0.0;

        _vis_acc.header.stamp = _odom.header.stamp;

        _vis_acc.points.clear();
        pt.x = _cmd.position.x;
        pt.y = _cmd.position.y;
        pt.z = _cmd.position.z;

        _vis_acc.points.push_back(pt);
        
        pt.x = _cmd.position.x + _cmd.acceleration.x;
        pt.y = _cmd.position.y + _cmd.acceleration.y;
        pt.z = _cmd.position.z + _cmd.acceleration.z;

        _vis_acc.points.push_back(pt);

        _vis_acc.scale.x = 0.1;
        _vis_acc.scale.y = 0.2;
        _vis_acc.scale.z = 0.2;

        _vis_acc_pub.publish(_vis_acc);
    }
};

int main(int argc, char ** argv)
{
    ros::init(argc, argv, "gradient_trajectory_server_node");
    ros::NodeHandle handle("~");

    handle.param("optimization/poly_order_min", _poly_order_min,  5);
    handle.param("optimization/poly_order_max", _poly_order_max,  10);
    TrajectoryServer server(handle);

    Bernstein _bernstein;
    if(_bernstein.setParam(_poly_order_min, _poly_order_max, 3) == -1) 
    {
        ROS_ERROR(" The trajectory order is set beyond the library's scope, please re-set ");
    }

    server.CList  = _bernstein.getC();
    server.CvList = _bernstein.getC_v();
    server.CaList = _bernstein.getC_a();

    sleep(1);
    ros::spin();

    return 0;
}

