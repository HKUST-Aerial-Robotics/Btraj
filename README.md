# Btraj
## 1.Introduction

Btraj is an online UAV plannning framework used to generate safe, dynamically feasible trajectories in previous unknown environemnts.

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
  Our testing environment: **Ubuntu** 16.04, **ROS** Kinetic.


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
  We use **mosek** for solving quadratic program(QP) and [sdf_tools](https://github.com/UM-ARM-Lab/sdf_tools) for building euclidean distance field.

## 7.Licence
The source code is released under [GPLv3](http://www.gnu.org/licenses/) license.
