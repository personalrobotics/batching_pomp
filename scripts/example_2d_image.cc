#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>

#include <boost/shared_ptr.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/function.hpp>
#include <boost/program_options.hpp>

#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/base/SpaceInformation.h>
#include <ompl/base/ProblemDefinition.h>
#include <ompl/base/Planner.h>
#include <ompl/geometric/PathGeometric.h>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <batching_pomp/BatchingPOMP.hpp>

namespace po = boost::program_options;
using batching_pomp::BatchingPOMP;

/// Check if the point is a black pixel or not
/// This is bound to the stateValidityChecker of the ompl StateSpace
/// \param[in] img The b/w image which represents the 2D map
/// \param[in] state The ompl state to check for validity
/// \return True if the pixel corresponding to the state is free, False otherwise
bool isPointValid(cv::Mat img, const ompl::base::State * state)
{
    double* values = state->as<
      ompl::base::RealVectorStateSpace::StateType>()->values;
    int r = img.rows, c = img.cols;

    cv::Point p((int)(values[0]*c),(int)(values[1]*r));
    
    if(img.at<uchar>(p.y,p.x)==0) //Black pixel
        return false;

    return true;
}

/// Create a 2D ompl scoped state from 2 floating point numbers
/// \param[in] space The ompl StateSpace for the problem
/// \param[in] x,y The floating point values of the two DOFs
/// \return The ompl RealVectorStateSpace (of dimension 2) that contains (x,y)
ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
make_state(const ompl::base::StateSpacePtr space, double x, double y)
{
   ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
      state(space);
   double * values = state->as<
      ompl::base::RealVectorStateSpace::StateType>()->values;
   values[0] = x;
   values[1] = y;
   return state;
}

/// Overlay the path on the map image
/// \param[in] img The OpenCV Mat structure containing the map image
/// \param[in] state_list The sequence of states along path
void displayPath(cv::Mat& img, std::vector< std::vector<double> > state_list)
{
    std::size_t nstates = state_list.size();
    cv::Mat copy = img.clone();
    int r=img.rows,c=img.cols;

    for(std::size_t i=0; i<nstates-1; i++){
        
        cv::Point p((int)((state_list[i][0])*c),(int)((state_list[i][1])*r));
        cv::Point p_next((int)((state_list[i+1][0])*c),(int)((state_list[i+1][1])*r));
        
        line(copy,p,p_next,cv::Scalar(255,0,0),3);
        circle(copy,p,6,cv::Scalar(0,0,255),-1);
    }
    cv::imshow("Path",copy);
    cv::waitKey(0);
}

/// Get ompl ScopedState from geometric path
/// \param[in] path The ompl geometric path to get states from
/// \param[in] idx The index of the state along the path
/// \return The ompl scoped state at path[idx]
ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
get_path_state(ompl::geometric::PathGeometric* path, size_t idx)
{
   ompl::base::StateSpacePtr space = path->getSpaceInformation()->getStateSpace();
   ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
      traj_state(space, path->getState(idx));
   return traj_state;
}

int main(int argc, char* argv[])
{

  po::options_description desc("2d Map Planner Options");

  desc.add_options()
      ("roadmapfile,f",po::value< std::string >()->required(), "Path to graph file")
      ("mapfile,m",po::value< std::string >()->required(),"Environment map image")
      ("start_coords,s",po::value<std::vector<float> >()->multitoken(),"start coordinates")
      ("goal_coords,g",po::value<std::vector<float> >()->multitoken(),"goal coordinates")
      ("vertices_type,v",po::value< std::string >()->required(),"Type of sequence to generate vertices (halton/random)")
      ("batching_type,b",po::value< std::string >()->required(), "Type of roadmap batching")
  ;

  // Read arguments
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string map_file(vm["mapfile"].as<std::string>());
  std::string graph_file(vm["roadmapfile"].as<std::string>());
  std::vector<float> start(vm["start_coords"].as<std::vector< float> >());
  std::vector<float> goal(vm["goal_coords"].as<std::vector< float> >());
  std::string vtype(vm["vertices_type"].as<std::string>());
  std::string btype(vm["batching_type"].as<std::string>());

  // Read image
  cv::Mat img = cv::imread(map_file,0); //BnW for computing
  cv::Mat img_disp = cv::imread(map_file,1); //Colour for display

  // 2-dimensional state space
  ompl::base::StateSpacePtr 
      space(new ompl::base::RealVectorStateSpace(2));
  space->as<ompl::base::RealVectorStateSpace>()->setBounds(0.0,1.0);
  space->setLongestValidSegmentFraction(
      0.01 / space->getMaximumExtent());
  space->setup();

  //Space info
  std::function<bool(const ompl::base::State*)> isStateValid = std::bind(
                          isPointValid,img,std::placeholders::_1);
  ompl::base::SpaceInformationPtr si(
      new ompl::base::SpaceInformation(space));
  si->setStateValidityChecker(isStateValid);
  si->setup();

  //Problem Definition
  ompl::base::ProblemDefinitionPtr pdef(
    new ompl::base::ProblemDefinition(si));
  pdef->addStartState(make_state(space, start[0], start[1]));
  pdef->setGoalState(make_state(space, goal[0], goal[1]));

  std::vector< std::vector<double> > state_list;

  // Set planner params
  // Alpha set to 0.2 and prune threshold set to 0.05 by default
  // Other params obtained from user options
  batching_pomp::BatchingPOMP bpPlanner(si);
  
  ompl::base::ParamSet& param_set = bpPlanner.params();
  param_set.setParam("increment","0.2");
  param_set.setParam("graph_type",vtype);
  param_set.setParam("batching_type",btype);
  param_set.setParam("prune_threshold","0.05");
  param_set.setParam("roadmap_filename",graph_file);

  bpPlanner.setup();
  bpPlanner.setProblemDefinition(pdef);

  ompl::base::PlannerStatus status;

  do
  {
    status = bpPlanner.solve(ompl::base::plannerNonTerminatingCondition());

    if(status == ompl::base::PlannerStatus::EXACT_SOLUTION 
      || status == ompl::base::PlannerStatus::APPROXIMATE_SOLUTION) {

      auto path = std::dynamic_pointer_cast<ompl::geometric::PathGeometric>(
            pdef->getSolutionPath());

      std::cout<<path->length()<<std::endl;

      std::size_t pathSize{path -> as< ompl::geometric::PathGeometric >()->getStateCount()};
      state_list.resize(pathSize);

      for(size_t i = 0; i < pathSize; i++) {
        ompl::base::ScopedState<ompl::base::RealVectorStateSpace> st = get_path_state(path -> as< ompl::geometric::PathGeometric >(),i);
        state_list[i] = st.reals();
      }

      displayPath(img_disp,state_list);
      //std::cin.ignore();
      pdef->clearSolutionPaths();
    }

  }while(status != ompl::base::PlannerStatus::TIMEOUT &&
    status != ompl::base::PlannerStatus::EXACT_SOLUTION);

  return 0;
}