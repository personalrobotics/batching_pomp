cmake_minimum_required(VERSION 2.8.12)
project(batching_pomp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#================================================================================
# Dependencies
find_package(catkin REQUIRED)
find_package(OMPL REQUIRED)
find_package(Boost REQUIRED COMPONENTS graph program_options)
find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)

catkin_package(
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
    DEPENDS Boost EIGEN3 OMPL
)

include_directories(
    include
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIRS}
    ${OMPL_INCLUDE_DIRS}
    ${catkin_INCLUDE_DIRS}
    ${OpenCV_INCLUDE_DIRS}
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

# For example script
add_executable(example_2d scripts/example_2d_image)
target_link_libraries(example_2d ${PROJECT_NAME} ${OpenCV_LIBS})

# For belief-resampling script
add_executable(example_2d_resample scripts/example_2d_resampling)
target_link_libraries(example_2d_resample ${PROJECT_NAME} ${OpenCV_LIBS})
