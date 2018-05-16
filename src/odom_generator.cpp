#include <iostream>
#include <math.h>
#include <random>
#include <eigen3/Eigen/Dense>
#include <ros/ros.h>
#include <ros/console.h>
#include <nav_msgs/Odometry.h>
#include <nav_msgs/Path.h>
#include <geometry_msgs/PoseStamped.h>
#include <visualization_msgs/MarkerArray.h>
#include <visualization_msgs/Marker.h>
#include <tf/tf.h>
#include <tf/transform_datatypes.h>
#include <tf/transform_broadcaster.h>
#include "quadrotor_msgs/PositionCommand.h"

ros::Subscriber _cmd_sub;
ros::Publisher  _odom_pub;

quadrotor_msgs::PositionCommand _cmd;
double _init_x, _init_y, _init_z;

bool rcv_cmd = false;
void rcvPosCmdCallBack(const quadrotor_msgs::PositionCommand cmd)
{	
	rcv_cmd = true;
	_cmd    = cmd;
}

void pubOdom()
{	
	nav_msgs::Odometry odom;
	odom.header.stamp    = ros::Time::now();
	odom.header.frame_id = "uav";

	if(rcv_cmd)
	{
	    odom.pose.pose.position.x = _cmd.position.x;
	    odom.pose.pose.position.y = _cmd.position.y;
	    odom.pose.pose.position.z = _cmd.position.z;

	    odom.twist.twist.linear.x = _cmd.velocity.x;
	    odom.twist.twist.linear.y = _cmd.velocity.y;
	    odom.twist.twist.linear.z = _cmd.velocity.z;

	    odom.twist.twist.angular.x = _cmd.acceleration.x;
	    odom.twist.twist.angular.y = _cmd.acceleration.y;
	    odom.twist.twist.angular.z = _cmd.acceleration.z;
	}
	else
	{
		odom.pose.pose.position.x = _init_x;
	    odom.pose.pose.position.y = _init_y;
	    odom.pose.pose.position.z = _init_z;

	    odom.twist.twist.linear.x = 0.0;
	    odom.twist.twist.linear.y = 0.0;
	    odom.twist.twist.linear.z = 0.0;

	    odom.twist.twist.angular.x = 0.0;
	    odom.twist.twist.angular.y = 0.0;
	    odom.twist.twist.angular.z = 0.0;
	}

    _odom_pub.publish(odom);
}

int main (int argc, char** argv) 
{        
    ros::init (argc, argv, "odom_generator");
    ros::NodeHandle nh( "~" );

    nh.param("init_x", _init_x,  0.0);
    nh.param("init_y", _init_y,  0.0);
    nh.param("init_z", _init_z,  0.0);

    _cmd_sub  = nh.subscribe( "command", 1, rcvPosCmdCallBack );
    _odom_pub = nh.advertise<nav_msgs::Odometry>("odometry", 1);                      

    ros::Rate rate(100);
    bool status = ros::ok();
    while(status) 
    {
		pubOdom();                   
        ros::spinOnce();
        status = ros::ok();
        rate.sleep();
    }

    return 0;
}