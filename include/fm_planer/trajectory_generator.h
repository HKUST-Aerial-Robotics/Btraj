#ifndef _TRAJECTORY_GENERATOR_H_
#define _TRAJECTORY_GENERATOR_H_

#include <stdio.h>
#include <ros/ros.h>
#include <ros/console.h>
#include <iostream>
#include <fstream>
#include <string>
#include "fm_planer/mosek.h"
#include "fm_planer/bezier_base.h"
#include "fm_planer/dataType.h"

using namespace std;
using namespace Eigen;

static void MSKAPI printstr(void *handle, MSKCONST char str[])
{
  printf("%s",str);
}

class TrajectoryGenerator {
private:

public:
        TrajectoryGenerator(){}
        ~TrajectoryGenerator(){}

        /* Use Bezier curve for the trajectory */
        MatrixXd BezierPloyCoeffGeneration(
            const vector<Cube> &corridor,
            const MatrixXd &MQM,
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
            const bool & isLimitAcc );  // define the order to which we minimize.   1 -- velocity, 2 -- acceleration, 3 -- jerk, 4 -- snap  

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