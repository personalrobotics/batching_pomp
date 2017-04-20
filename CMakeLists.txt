cmake_minimum_required(VERSION 2.8.12)
project(batching_pomp)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#================================================================================
# Dependencies
#

find_package(catkin REQUIRED COMPONENTS cmake_modules)
find_package(FLANN REQUIRED)
find_package(OMPL REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)
find_package(Eigen REQUIRED)

# TODO - Figure out if ROS Kinetic using std::shared_ptr affects

catkin_package(
   INCLUDE_DIRS include ${CATKIN_DEVEL_PREFIX}/include
   LIBRARIES ${PROJECT_NAME}
   DEPENDS Boost Eigen
)

# TODO - Set CMAKE_BUILD_TYPE to Release?

include_directories(
    include
    ${CATKIN_DEVEL_PREFIX}/include
    ${catkin_INCLUDE_DIRS}
    ${OMPL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${Eigen_INCLUDE_DIRS}
)

link_libraries(
   ${catkin_LIBRARIES}
   ${OMPL_LIBRARIES}
   ${Boost_LIBRARIES} )

add_definitions("-Wall")

#add_library(${PROJECT_NAME})
#    
#)

add_executable(model_tester ignore/test_model.cpp)
target_link_libraries(model_tester ${catkin_LIBRARIES} ${OMPL_LIBRARIES})
