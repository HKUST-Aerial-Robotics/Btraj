#include <iostream>
#include <ros/ros.h>
#include <ros/console.h>
#include <Eigen/Eigen>
#include <math.h>
#include "backward.hpp"
#include "data_type.h"
#include <sdf_tools/collision_map.hpp>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

class kinoGridPathFinder
{
	private:
		inline Eigen::Vector3d gridIndex2coord(Eigen::Vector3i index);
		inline Eigen::Vector3i coord2gridIndex(Eigen::Vector3d pt);
		inline void coord2gridIndexFast(double x, double y, double z, int & id_x, int & id_y, int & id_z);
		inline KinoGridNodePtr pos2KinoGridNodePtr(Eigen::Vector3d pos);
		inline double getHeu(KinoGridNodePtr node1, KinoGridNodePtr node2);

		std::vector<KinoGridNodePtr> retrievePath(KinoGridNodePtr current);

		double resolution, inv_resolution;
		double gl_xl, gl_yl, gl_zl;
		double gl_xu, gl_yu, gl_zu;
		double tie_breaker = 1.0 + 1.0 / 10000;

		double t_prop  = 0.5;
		double t_delta = 0.1;
		std::vector<double> t_steps;

		double u_max   = 1.0;
		double u_delta = 0.5;

		double w_h   = 1.0;
		double w_t   = 1.0;
		double v_max = 2.0;

		double time_in_forward;
		int num_ope = 0;

		std::vector<KinoGridNodePtr> expandedNodes;
		std::vector<KinoGridNodePtr> gridPath;
		KinoGridNodePtr terminatePtr;

		int GLX_SIZE, GLY_SIZE, GLZ_SIZE;
		int X_SIZE, Y_SIZE, Z_SIZE;

		KinoGridNodePtr *** KinoGridNodeMap;
		std::multimap<double, KinoGridNodePtr> openSet;

	public:
		kinoGridPathFinder( Eigen::Vector3i GL_size, Eigen::Vector3i LOC_size)
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
		kinoGridPathFinder(){};
		~kinoGridPathFinder(){};

		void setParameter(double t_ptop_, double t_delta_, double u_max_, double u_delta_, double w_h_, double w_t_, double v_max_ );
		void initGridNodeMap(double _resolution, Eigen::Vector3d global_xyz_l, Eigen::Vector3d global_xyz_u);
		void linkLocalMap(sdf_tools::CollisionMapGrid * local_map, Eigen::Vector3d xyz_l);
		void hybridAstarSearch(Eigen::Vector3d start_pt, Eigen::Vector3d start_vel, Eigen::Vector3d end_pt, Eigen::Vector3d end_vel);
		inline bool forwardSimulation(KinoGridNodePtr p_cur, Eigen::Vector3d u, KinoGridNodePtr succPtr);
		inline void getState( Eigen::VectorXd & xt, Eigen::Vector3d u, double t );
		std::vector<Eigen::Vector3d> getKinoTraj(double resolution);

		void resetLocalMap();
		void resetPath();

		std::vector<Eigen::Vector3d> getPath();
		std::vector<Eigen::Vector3d> getVisitedNodes();
};