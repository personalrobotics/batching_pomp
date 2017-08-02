cmake_minimum_required(VERSION 2.8.12)
project(batching_pomp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#================================================================================
# Dependencies
find_package(catkin REQUIRED)
find_package(OMPL REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph)
find_package(Eigen3 REQUIRED)

catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
    DEPENDS Boost Eigen OMPL
)

include_directories(
    include
    ${CATKIN_DEVEL_PREFIX}/include
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${OMPL_INCLUDE_DIRS}
    ${catkin_INCLUDE_DIRS}
)
link_directories(
    ${OMPL_LIBRARY_DIRS}
    ${catkin_LIBRARY_DIRS}
)
add_definitions(
    ${EIGEN3_DEFINITIONS}
)


add_library(${PROJECT_NAME}
    src/BatchingPOMP.cpp
)

target_link_libraries(${PROJECT_NAME}
    ${Boost_LIBRARIES}
    ${OMPL_LIBRARIES}
)

install(TARGETS batching_pomp
    LIBRARY DESTINATION "${CATKIN_PACKAGE_LIB_DESTINATION}"
)
install(DIRECTORY "include/${PROJECT_NAME}/"
    DESTINATION "${CATKIN_PACKAGE_INCLUDE_DESTINATION}"
    PATTERN "*.in" EXCLUDE
    PATTERN ".svn" EXCLUDE
)
install(DIRECTORY
    "${CATKIN_DEVEL_PREFIX}/include/${PROJECT_NAME}/"
    DESTINATION "${CATKIN_PACKAGE_INCLUDE_DESTINATION}"
)