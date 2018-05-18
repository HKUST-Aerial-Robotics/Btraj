# Btraj
## 1.Introduction

Btraj is an online UAV plannning framework used to generate safe, dynamically feasible trajectories in previous unknown environemnts. It can be divided as front-end path finding module and back-end trajectory optimization module. In the front-end, we provide two alternates: Fast Marching*(FM*) on a velocity field and A* on a pure grid map. A flight corridor conststs of cubes are generated based on the path. In the back-end, we utilize properties of Bezier curve to confine the piecewise Bezier curves entirely within the corridor and dynamical limits. 

**Authors:**[Fei Gao](https://ustfei.com/) and [Shaojie Shen](http://www.ece.ust.hk/ece.php/profile/facultydetail/eeshaojie) from the [HUKST Aerial Robotics Group](uav.ust.hk).

**Disclaimer**

This is research code, any fitness for a particular purpose is disclaimed.

**Related Paper**
* **Online Safe Trajectory Generation For Quadrotors
Using Fast Marching Method and Bernstein Basis Polynomial,** Fei Gao, William Wu, Yi Lin and Shaojie Shen

Video of this paper can be found:

<a href="https://www.youtube.com/watch?v=Dn6pXL3GqeY" target="_blank"><img src="https://img.youtube.com/vi/Dn6pXL3GqeY/0.jpg" 
alt="video" width="540" height="360" border="10" /></a>


If you use this planning framework for your academic research, please cite our related paper.
```
@inproceedings{Fei2018ICRA,
	Address = {Brisbane, Australia},
	Author = {F. Gao and W.Wu and Y. Lin and S. Shen},
	Booktitle = {Online Safe Trajectory Generation For Quadrotors
Using Fast Marching Method and Bernstein Basis Polynomial},
	Title = {Proc. of the {IEEE} Intl. Conf. on Robot. and Autom.},
	Month = May,
	Year = {2018}}
}
```
## 2.Prerequisities
- Our testing environment: **Ubuntu** 16.04, **ROS** Kinetic.
- We provide a simple simulation to test the code. To run the simulation, you should install [armadillo](http://arma.sourceforge.net/), which is a c++ linear algebra library. Then clone and compile [plan_utils](https://github.com/HKUST-Aerial-Robotics/plan_utils), which contains several ROS-package used for running the simulation.
```
  sudo apt-get install libarmadillo-dev
  cd ~/catkin_ws/src
  git clone https://github.com/HKUST-Aerial-Robotics/plan_utils.git
  cd ../
  catkin_make
  source ~/catkin_ws/devel/setup.bash
```

## 3.Build on ROS
  Clone the repository to your catkin workspace and catkin_make. For example:
```
  cd ~/catkin_ws/src
  git clone https://github.com/HKUST-Aerial-Robotics/Btraj.git
  cd ../
  catkin_make
  source ~/catkin_ws/devel/setup.bash
```

## 4.Install Mosek


## 4.Usage


## 6.Acknowledgements
  We use **mosek** for solving quadratic program(QP), [fast_methods](https://github.com/jvgomez/fast_methods) for performing general fast marching method and [sdf_tools](https://github.com/UM-ARM-Lab/sdf_tools) for building euclidean distance field.

## 7.Licence
The source code is released under [GPLv3](http://www.gnu.org/licenses/) license.

## 8.Notes
- The code has not been deeply tested, if you find any problems, do not hesitate to raise a issue or write e-mail to me directly.
- The code is written for research purpose and has not been fully optimized. In the future I will add more functionalities and improve efficiency, and also add more comment. 

