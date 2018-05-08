#include "trajectory_generator_waypoint.h"
#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;    
using namespace Eigen;

#define inf 1>>30

TrajectoryGeneratorWaypoint::TrajectoryGeneratorWaypoint(){}

TrajectoryGeneratorWaypoint::~TrajectoryGeneratorWaypoint(){}

Eigen::MatrixXd TrajectoryGeneratorWaypoint::PolyQPGeneration(
            const Eigen::MatrixXd &Path,
            const Eigen::Vector3d &Vel,
            const Eigen::Vector3d &Acc,
            const Eigen::VectorXd &Time) 
{     
      /*   Get initial trajectory which is a straight line( 0 zero end velocity and acceleration ) minimum snap trajectory or truly minumum snap trajectory  */
      /*ros::Time time_1 = ros::Time::now();
      ros::Time time_2, time_3;*/

      int m = Time.size();
      MatrixXd PolyCoeff(m, 3 * 6);
      VectorXd Px(6 * m), Py(6 * m), Pz(6 * m);

      int num_f, num_p; // number of fixed and free variables
      int num_d;        // number of all segments' derivatives
      const static auto Factorial = [](int x){
          int fac = 1;

          for(int i = x; i > 0; i--)
              fac = fac * i;
            
          return fac;
      };

      /*   Produce Mapping Matrix A to the entire trajectory.   */
      MatrixXd Ab;
      MatrixXd A = MatrixXd::Zero(m * 6, m * 6);

      for(int k = 0; k < m; k++){
          Ab = Eigen::MatrixXd::Zero(6, 6);
          for(int i = 0; i < 3; i++){
              Ab(2 * i, i) = Factorial(i);
              for(int j = i; j < 6; j++)
                  Ab( 2 * i + 1, j ) = (double)Factorial(j) / (double)Factorial( j - i ) * pow( Time(k), j - i );
          }

          A.block(k * 6, k * 6, 6, 6) = Ab;    
      }
      
      MatrixXd A_inv   = A.inverse();
      //ROS_WARN("[Generator] A finished");
      /*   Produce the dereivatives in X, Y and Z axis directly.  */
      VectorXd Dx = VectorXd::Zero(m * 6);
      VectorXd Dy = VectorXd::Zero(m * 6);
      VectorXd Dz = VectorXd::Zero(m * 6);

      for(int k = 1; k < (m + 1); k ++ ){
          Dx((k-1)*6) = Path(k - 1, 0); Dx((k-1)*6 + 1) = Path(k, 0); 
          Dy((k-1)*6) = Path(k - 1, 1); Dy((k-1)*6 + 1) = Path(k, 1); 
          Dz((k-1)*6) = Path(k - 1, 2); Dz((k-1)*6 + 1) = Path(k, 2); 
          
          if( k == 1){
              Dx((k-1)*6 + 2) = Vel(0);
              Dy((k-1)*6 + 2) = Vel(1); 
              Dz((k-1)*6 + 2) = Vel(2);

              Dx((k-1)*6 + 4) = Acc(0);
              Dy((k-1)*6 + 4) = Acc(1); 
              Dz((k-1)*6 + 4) = Acc(2);
          }
      }

      /*   Produce the Minimum Snap cost function, the Hessian Matrix   */
      MatrixXd H = MatrixXd::Zero( m * 6, m * 6 );
      
      for(int k = 0; k < m; k ++){
          for(int i = 3; i < 6; i ++){
              for(int j = 3; j < 6; j ++){
                  H( k*6 + i, k*6 + j ) = (double)i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) / (double)(i + j - 5) * pow( Time(k), (i + j - 5) );
              }
          }
      }

/*      MatrixXd H1d = MatrixXd::Zero( 6, 6 );
      cout<<"Time(0): "<<Time(0)<<endl;
      for(int i = 3; i < 6; i ++){
          for(int j = 3; j < 6; j ++){
              H1d( i, j ) = (double)i * (i - 1) * (i - 2) * j * (j - 1) * (j - 2) * pow( Time(0), (i + j - 5) ) / (double)(i + j - 5);
          }
      }
      cout<<"H1d:\n"<<H1d<<endl;*/

      _Q = H; // Now only minumum snap is used in the cost

      if( m > 1)
      {   
          //ROS_WARN("[Generator] H finished");
          MatrixXd Ct; // The transpose of selection matrix C
          MatrixXd C;  // The selection matrix C

          num_f = 2 * m + 4; //3 + 3 + (m - 1) * 2 = 2m + 4
          num_p = 2 * m - 2; //(m - 1) * 2 = 2m - 2
          num_d = 6 * m;
          Ct = MatrixXd::Zero(num_d, num_f + num_p); 
          Ct( 0, 0 ) = 1; Ct( 2, 1 ) = 1;         Ct( 4, 2 ) = 1; // stack the start point
          Ct( 1, 3 ) = 1; Ct( 3, 2 * m + 4 ) = 1; Ct( 5, 2 * m + 5 ) = 1; 

          Ct(6 * (m - 1) + 0, 2 * m + 0) = 1; 
          Ct(6 * (m - 1) + 1, 2 * m + 1) = 1; // Stack the end point
          Ct(6 * (m - 1) + 2, 4 * m + 0) = 1;
          Ct(6 * (m - 1) + 3, 2 * m + 2) = 1; // Stack the end point
          Ct(6 * (m - 1) + 4, 4 * m + 1) = 1;
          Ct(6 * (m - 1) + 5, 2 * m + 3) = 1; // Stack the end point

          for(int j = 2; j < m; j ++ ){
              Ct( 6 * (j - 1) + 0, 2 + 2 * (j - 1)         + 0 ) = 1;
              Ct( 6 * (j - 1) + 1, 2 + 2 * (j - 1)         + 1 ) = 1;
              Ct( 6 * (j - 1) + 2, 2 * m + 4 + 2 * (j - 2) + 0 ) = 1;
              Ct( 6 * (j - 1) + 3, 2 * m + 4 + 2 * (j - 1) + 0 ) = 1;
              Ct( 6 * (j - 1) + 4, 2 * m + 4 + 2 * (j - 2) + 1 ) = 1;
              Ct( 6 * (j - 1) + 5, 2 * m + 4 + 2 * (j - 1) + 1 ) = 1;
          }

          C = Ct.transpose();
          MatrixXd A_invC  = A_inv * Ct;

          //ROS_WARN("[Generator] in type 1, C finished");
          VectorXd Dx1 = C * Dx;
          VectorXd Dy1 = C * Dy;
          VectorXd Dz1 = C * Dz;

          /*time_2 = ros::Time::now();
          ROS_WARN("[Waypoint-Traj Solver] time in stacking variables is %f", (time_2 - time_1).toSec() );*/
          //ROS_WARN("[Generator] case segment > 1");
          MatrixXd R   = A_invC.transpose() * _Q *  A_invC;
          
          //ROS_WARN("[Generator] R finished");
          VectorXd Dxf(2 * m + 4), Dyf(2 * m + 4), Dzf(2 * m + 4);
          
          Dxf = Dx1.segment( 0, 2 * m + 4 );
          Dyf = Dy1.segment( 0, 2 * m + 4 );
          Dzf = Dz1.segment( 0, 2 * m + 4 );

          MatrixXd Rff(2 * m + 4, 2 * m + 4);
          MatrixXd Rfp(2 * m + 4, 2 * m - 2);
          MatrixXd Rpf(2 * m - 2, 2 * m + 4);
          MatrixXd Rpp(2 * m - 2, 2 * m - 2);
     
          Rff = R.block(0, 0, 2 * m + 4, 2 * m + 4);
          Rfp = R.block(0, 2 * m + 4, 2 * m + 4, 2 * m - 2);
          Rpf = R.block(2 * m + 4, 0,         2 * m - 2, 2 * m + 4);
          Rpp = R.block(2 * m + 4, 2 * m + 4, 2 * m - 2, 2 * m - 2);

          MatrixXd Rpp_inv = Rpp.inverse();
          //ROS_WARN("[Generator] R blocks finished");
          VectorXd Dxp(2 * m - 2), Dyp(2 * m - 2), Dzp(2 * m - 2);
          Dxp = - (Rpp_inv * Rfp.transpose()) * Dxf;
          Dyp = - (Rpp_inv * Rfp.transpose()) * Dyf;
          Dzp = - (Rpp_inv * Rfp.transpose()) * Dzf;

          //ROS_WARN("[Generator] Dp blocks finished");
          Dx1.segment(2 * m + 4, 2 * m - 2) = Dxp;
          Dy1.segment(2 * m + 4, 2 * m - 2) = Dyp;
          Dz1.segment(2 * m + 4, 2 * m - 2) = Dzp;

          //ROS_WARN("[Generator] Dx1 resembled finished");
          Px = A_invC * Dx1;
          Py = A_invC * Dy1;
          Pz = A_invC * Dz1;

          //cout<<"J: "<<Dx1.transpose() * R * Dx1 + Dy1.transpose() * R * Dy1 + Dz1.transpose() * R * Dz1<<endl;
      }
      else
      {   
          //ROS_WARN("[Generator] case segment = 1");
          Px = A_inv * Dx;
          Py = A_inv * Dy;
          Pz = A_inv * Dz;
          //time_2 = ros::Time::now();
      }

      _Px = Px;
      _Py = Py;
      _Pz = Pz;

      for(int i = 0; i < m; i ++)
      {
          PolyCoeff.block(i, 0,  1, 6) = Px.segment( i * 6, 6 ).transpose();
          PolyCoeff.block(i, 6,  1, 6) = Py.segment( i * 6, 6 ).transpose();
          PolyCoeff.block(i, 12, 1, 6) = Pz.segment( i * 6, 6 ).transpose();
      }

      /*time_3 = ros::Time::now();
      ROS_WARN("[Waypoint-Traj Solver] time in doing calculation is %f", (time_3 - time_2).toSec() );*/
      //ROS_WARN("[Generator] Unconstrained QP solved");

      //cout<<"objective: "<<_Px.transpose() * _Q * _Px + _Py.transpose() * _Q * _Py + _Pz.transpose() * _Q * _Pz<<endl;
      return PolyCoeff;
}  

double TrajectoryGeneratorWaypoint::getObjective()
{ 
      _qp_cost = (_Px.transpose() * _Q * _Px + _Py.transpose() * _Q * _Py + _Pz.transpose() * _Q * _Pz)(0);
      return _qp_cost; 
}
