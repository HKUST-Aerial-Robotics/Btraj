#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include "backward.hpp"
#include "gridNode.h"

using namespace std;
using namespace Eigen;

class gridPathFinder
{
	private:
		Vector3d    gridIndex2coord(Vector3i index);
		Vector3i    coord2gridIndex(Vector3d pt);
		GridNodePtr pos2gridNodePtr(Vector3d pos);
		double getDiagHeu(GridNodePtr node1, GridNodePtr node2);
		double getManhHeu(GridNodePtr node1, GridNodePtr node2);
		double getEuclHeu(GridNodePtr node1, GridNodePtr node2);
		double getHeu(GridNodePtr node1, GridNodePtr node2);

		vector<GridNodePtr> retrievePath(GridNodePtr current);

		double resolution;
		double gl_xl, gl_yl, gl_zl;
		double tie_breaker = 1.0 + 1.0 / 10000;

		vector<GridNodePtr> expandedNodes;
		vector<GridNodePtr> gridPath;

		int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
		int X_SIZE, Y_SIZE, Z_SIZE;
		Vector3d loc_map_o;

		int *** gridMap;
		GridNodePtr *** GridNodeMap;

		multimap<double, GridNodePtr> openSet;
	public:
		gridPathFinder( Vector3i GL_size, Vector3i LOC_size)
		{	
			// size of a big big global grid map
			GLX_SIZE = GL_size(0);
			GLY_SIZE = GL_size(1);
			GLZ_SIZE = GL_size(2);

			// size of local map, for local obs recording
			X_SIZE = LOC_size(0);
			Y_SIZE = LOC_size(1);
			Z_SIZE = LOC_size(2);
		};
		~gridPathFinder(){};

		void initGridNodeMap(double _resolution, Vector3d global_xyz_l);
		void linkLocalMap(int *** map, Vector3d xyz_l);
		void resetLocalMap();
		bool CheckGuidePathCollision();
		void AstarSearch(Vector3d start_pt, Vector3d end_pt);
		void resetPath();

		vector<GridNodePtr> getPath();
		vector<GridNodePtr> getVisitedNodes();
};