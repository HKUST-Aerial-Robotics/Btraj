#ifndef _DATA_TYPE_
#define _DATA_TYPE_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <vector>

#define _inf 1>>30
#define _PI M_PI

using namespace std;

struct GridNode;

typedef GridNode* GridNodePtr;

struct Idx
{
      int idx, idy, idz;

      Idx(int id_x, int id_y, int id_z)
      {
         idx = id_x;
         idy = id_y;
         idz = id_z;
      }

      Idx(){};
      ~Idx(){};
};
  
struct GridNode
{     
   int id; // 1--> open set, -1 --> closed set
   Eigen::Vector3d real_coord;
   Idx grid_index;
   double gScore, fScore;
   double occupied; //consistent with occupancy grid
   GridNodePtr cameFrom;
   std::multimap<double, GridNodePtr>::iterator nodeMapIt;

   GridNode(Idx index, Eigen::Vector3d coord)
   {  
      id = 0;
      grid_index = index;
      real_coord = coord;
      gScore = _inf;
      fScore = _inf;
      cameFrom = NULL;
   }

   GridNode(Idx index)
   {
      id = 0;
      grid_index = index;
      gScore = _inf;
      fScore = _inf;
      cameFrom = NULL;
   }

   GridNode(){};
   
   ~GridNode(){};
};

struct Cube;
struct Cube
{     
      Eigen::Vector3d p1, p2, p3, p4, p5, p6, p7, p8;   // the 8 vertex of a cube 
      Eigen::Vector3d pc; // the center of the cube
      bool valid;    // indicates whether this cube should be deleted

      double t; // time allocated to this cube
      vector< pair<double, double> > box;

      Eigen::Vector3d start_pt, end_pt; // the start point and end point of the object's trajectory in this cube
      Eigen::MatrixXd ctrlPts; // the object's trajectory which the drone should track in this cube
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

      // create a cube using 8 vertex and the center point
      Cube( Eigen::Vector3d p1_, 
            Eigen::Vector3d p2_, 
            Eigen::Vector3d p3_, 
            Eigen::Vector3d p4_, 
            Eigen::Vector3d p5_, 
            Eigen::Vector3d p6_, 
            Eigen::Vector3d p7_, 
            Eigen::Vector3d p8_,
            Eigen::Vector3d pc_)
      {
            p1 = p1_; p2 = p2_; p3 = p3_; p4 = p4_;
            p5 = p5_; p6 = p6_; p7 = p7_; p8 = p8_;
            pc = pc_;
            valid = true;
      
            t = 0.0;
            box.resize(3);
            ctrlPts = Eigen::MatrixXd::Zero(3, 3);
      }

      // create a inscribe cube of a ball using the center point and the radius of the ball
      void setVertex(vector<double> pos_1, vector<double> pos_2, vector<double> pos_3, vector<double> pos_4,
                     vector<double> pos_5, vector<double> pos_6, vector<double> pos_7, vector<double> pos_8, double res )
      {
            p1 = Eigen::Vector3d(pos_1[0] + res, pos_1[1] - res, pos_1[2] + res);  
            p2 = Eigen::Vector3d(pos_2[0] + res, pos_2[1] + res, pos_2[2] + res);  
            p3 = Eigen::Vector3d(pos_3[0] - res, pos_3[1] + res, pos_3[2] + res);  
            p4 = Eigen::Vector3d(pos_4[0] - res, pos_4[1] - res, pos_4[2] + res);  

            p5 = Eigen::Vector3d(pos_5[0] + res, pos_5[1] - res, pos_5[2] - res);  
            p6 = Eigen::Vector3d(pos_6[0] + res, pos_6[1] + res, pos_6[2] - res);  
            p7 = Eigen::Vector3d(pos_7[0] - res, pos_7[1] + res, pos_7[2] - res);  
            p8 = Eigen::Vector3d(pos_8[0] - res, pos_8[1] - res, pos_8[2] - res);  

            /*pc(0) = (p1(0) + p4(0)) / 2.0; 
            pc(1) = (p1(1) + p2(1)) / 2.0; 
            pc(2) = (p1(2) + p5(2)) / 2.0; */
            setBox();
      }
      
      void setBox()
      {
            box.clear();
            box.resize(3);
            box[0] = make_pair( p4(0), p1(0) );
            box[1] = make_pair( p1(1), p2(1) );
            box[2] = make_pair( p5(2), p1(2) );
      }

      void printBox()
      {
            cout<<"pc:  "<<pc.transpose()<<endl;
            cout<<"p1:  "<<p1.transpose()<<endl;
            cout<<"p2:  "<<p2.transpose()<<endl;
            cout<<"p3:  "<<p3.transpose()<<endl;
            cout<<"p4:  "<<p4.transpose()<<endl;
            cout<<"p5:  "<<p5.transpose()<<endl;
            cout<<"p6:  "<<p6.transpose()<<endl;
            cout<<"p7:  "<<p7.transpose()<<endl;
            cout<<"p8:  "<<p8.transpose()<<endl;
      }

      Cube()
      {  
         pc = p1 = p2 = p3 = p4 = p5 = p6 = p7 = p8 = Eigen::VectorXd::Zero(3);
         valid = true;
         t = 0.0;
         box.resize(3);
         ctrlPts = Eigen::MatrixXd::Zero(3, 3);
      }

      ~Cube(){}
};

#endif