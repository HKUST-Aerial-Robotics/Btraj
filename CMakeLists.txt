cmake_minimum_required(VERSION 2.8.3)
project(bezier_planer)

find_package(catkin REQUIRED COMPONENTS rospy std_msgs sensor_msgs visualization_msgs tf image_transport cv_bridge arc_utilities message_generation)
find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)

set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
find_package(OpenCV REQUIRED)

catkin_package(
    INCLUDE_DIRS 
    include 
    LIBRARIES
    CATKIN_DEPENDS 
    rospy 
    std_msgs 
    sensor_msgs 
    visualization_msgs 
    image_transport 
    cv_bridge 
    arc_utilities 
    message_runtime 
    DEPENDS 
    Eigen3 
    OpenCV
    PCL)

link_directories( /home/bbgf/catkin_ws/src/bezier_planer/lib/mosek8 )

include_directories(
    include 
    SYSTEM 
    /home/bbgf/catkin_ws/src/bezier_planer/third_party/
    /home/bbgf/catkin_ws/devel/include/
    /home/bbgf/catkin_ws/src/sdf_tools/include/
    ${catkin_INCLUDE_DIRS} 
    ${Eigen3_INCLUDE_DIRS} 
    ${PCL_INCLUDE_DIRS}
    )

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -O3 -Wall") # -Wextra -Werror

# Simple test node for voxel_grid
add_executable(
    fm_trajectory_node 
    src/fm_trajectory_node.cpp
    src/bezier_base.cpp
    src/trajectory_generator.cpp
    src/trajectory_generator_socp.cpp
    src/trajectory_generator_waypoint.cpp
    third_party/fast_methods/console/console.cpp
    third_party/fast_methods/ndgridmap/cell.cpp
    third_party/fast_methods/fm/fmdata/fmcell.cpp
)
target_link_libraries(fm_trajectory_node ${catkin_LIBRARIES} sdf_tools ${PCL_LIBRARIES} mosek64)

add_executable ( fm_traj_server src/traj_server.cpp src/bezier_base.cpp)
target_link_libraries( fm_traj_server
                        ${catkin_LIBRARIES}
)

add_executable ( random_forest_sensing src/random_forest_sensing.cpp )
target_link_libraries( random_forest_sensing
                        ${catkin_LIBRARIES}
                        ${PCL_LIBRARIES}
)