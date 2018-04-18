/*

This header file is used to provide some basic mathematic support for the Bernstein-basis trajectory generation optimization problem. Includes:
1-: Mapping matrix maps the coefficients of the Bernstein basis (ie. control points) to Monomial basis. The mapping matrix range from order 3 to order 10
2-: Modulus list of the Bernstein basis to a given order. That is, pre-compute the constant-modulus (the 'n choose k' combinatorial) of the basis vector. 
	To save computation cost of frequently call this value. 

The class should be initialized to a instance before the trajectory generator called. 
Several initializer are provided, and the instance is initialized according to the given order of the control points.

*/
#ifndef _BEZIER_BASE_H_
#define _BEZIER_BASE_H_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace std;
using namespace Eigen;

class Bernstein
{
	private:
		
		vector<MatrixXd> S_foreList, S_backList;
		vector<MatrixXd> MQMList, MList;
		vector<VectorXd> CList, CvList, CaList, CjList;

		int _order_min, _order_max;  // The order of the polynomial in each segment, also the number of control points used in each segment
		int _min_order;              // The order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap    

	public:
		Bernstein(); // Empty constructor

		Bernstein(int poly_order_min, int poly_order_max, int min_order);

		~Bernstein();

		int setParam(int poly_order_min, int poly_order_max, int min_order);

		MatrixXd getS_fore(int order, double t);
		MatrixXd getS_back(int order, double t);

		vector<MatrixXd> getM();
		vector<MatrixXd> getMQM();	
		vector<VectorXd> getC();
		vector<VectorXd> getC_v();
		vector<VectorXd> getC_a();
		vector<VectorXd> getC_j();
};

#endif