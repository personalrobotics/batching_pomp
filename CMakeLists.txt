cmake_minimum_required(VERSION 2.8.12)
project(batching_pomp)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#================================================================================
# Dependencies


find_package(catkin REQUIRED COMPONENTS cmake_modules)
find_package(FLANN REQUIRED)
find_package(OMPL REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph program_options)
find_package(Eigen3 REQUIRED)

# TODO - Figure out if ROS Kinetic using std::shared_ptr affects

catkin_package(
   INCLUDE_DIRS include ${CATKIN_DEVEL_PREFIX}/include
   LIBRARIES ${PROJECT_NAME}
   DEPENDS Boost Eigen
)

set(CMAKE_BUILD_TYPE Release)

include_directories(
    include
    ${CATKIN_DEVEL_PREFIX}/include
    ${catkin_INCLUDE_DIRS}
    ${OMPL_INCLUDE_DIRS}
    ${Boost_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
)

link_libraries(
   ${catkin_LIBRARIES}
   ${OMPL_LIBRARIES}
   ${Boost_LIBRARIES} )

add_library(${PROJECT_NAME}
	src/BatchingPOMP.cpp)


add_executable(test_model ignore/test_model.cpp)
target_link_libraries(test_model ${PROJECT_NAME} ${OMPL_LIBRARIES})

add_executable(bare_planner ignore/test_bare_planner.cpp)
target_link_libraries(bare_planner ${PROJECT_NAME} ${OMPL_LIBRARIES} ${Boost_LIBRARIES} ${OpenCV_LIBS})

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)