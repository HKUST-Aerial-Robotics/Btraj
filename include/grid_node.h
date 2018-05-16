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

#define inf 1>>30

struct GridNode;

typedef GridNode* GridNodePtr;

struct GridNode
{     
   int id;        // 1--> open set, -1 --> closed set
   Eigen::Vector3d coord;
   Eigen::Vector3i index;
   
   double gScore, fScore;
   GridNodePtr cameFrom;
   multimap<double, GridNodePtr>::iterator nodeMapIt;
   double occupancy; 

   std::vector<GridNodePtr> hisNodeList; // use a list to record nodes in its history

   GridNode(Eigen::Vector3i _index)
   {  
      id = 0;
      index = _index;
      
      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
   }

   GridNode(Eigen::Vector3i _index, Eigen::Vector3d _coord)
   {  
      id = 0;
      index = _index;
      coord = _coord;

      gScore = inf;
      fScore = inf;
      cameFrom = NULL;
   }

   GridNode(){};
   
   ~GridNode(){};
};

#endif