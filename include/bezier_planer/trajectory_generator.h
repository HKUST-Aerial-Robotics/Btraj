#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "bezier_planer/mosek.h"
#include "bezier_planer/bezier_base.h"
#include "bezier_planer/dataType.h"

using namespace std;
using namespace Eigen;

class TrajectoryGenerator {
private:

public:
        TrajectoryGenerator(){}
        ~TrajectoryGenerator(){}

        /* Use Bezier curve for the trajectory */
       int BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const int minimize_order,
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc,
            double & obj,
            MatrixXd & PolyCoeff); 

        MatrixXd BezierPloyCoeffGenerationSOCP(
            const vector<Cube> &corridor,
            const MatrixXd &FM,
            const MatrixXd &pos,
            const MatrixXd &vel,
            const MatrixXd &acc,
            const double maxVel,
            const double maxAcc,
            const int traj_order,
            const int minimize_order,
            double & obj, 
            const double margin,
            const bool & isLimitVel,
            const bool & isLimitAcc );
};

#endif