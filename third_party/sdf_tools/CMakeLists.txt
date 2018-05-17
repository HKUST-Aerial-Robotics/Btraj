cmake_minimum_required(VERSION 2.8.3)
project(sdf_tools)

find_package(catkin REQUIRED COMPONENTS rospy std_msgs sensor_msgs visualization_msgs image_transport cv_bridge message_generation)
find_package(cmake_modules REQUIRED)
find_package(Eigen3 REQUIRED)
set(Eigen3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR})
find_package(OpenCV REQUIRED)

add_message_files(DIRECTORY msg FILES SDF.msg CollisionMap.msg)
generate_messages(DEPENDENCIES geometry_msgs std_msgs)

catkin_package(
    INCLUDE_DIRS include 
    LIBRARIES ${PROJECT_NAME} 
    CATKIN_DEPENDS 
    rospy 
    std_msgs 
    sensor_msgs 
    visualization_msgs 
    image_transport 
    cv_bridge 
    message_runtime 
    DEPENDS 
    Eigen3 
    OpenCV)

include_directories(include SYSTEM ${catkin_INCLUDE_DIRS} ${Eigen3_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS} -flto -O3 -Wall -Wextra -Werror")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

# SDF library
add_library(${PROJECT_NAME}
    include/${PROJECT_NAME}/collision_map.hpp
    include/${PROJECT_NAME}/dynamic_spatial_hashed_collision_map.hpp
    include/${PROJECT_NAME}/sdf.hpp
    src/${PROJECT_NAME}/collision_map.cpp
    src/${PROJECT_NAME}/dynamic_spatial_hashed_collision_map.cpp
    src/${PROJECT_NAME}/sdf.cpp)
add_dependencies(${PROJECT_NAME} ${catkin_EXPORTED_TARGETS} ${PROJECT_NAME}_gencpp)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})