#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <functional>
#include <chrono>
#include <random>

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

#include <batching_pomp/belief_resampling/SGDBasedResampler.hpp>
#include <batching_pomp/BatchingPOMP.hpp>
#include <batching_pomp/GraphTypes.hpp>

namespace po = boost::program_options;
using batching_pomp::BatchingPOMP;

/// Check if the point is a black pixel or not
/// This is bound to the stateValidityChecker of the ompl StateSpace
/// \param[in] img The b/w image which represents the 2D map
/// \param[in] state The ompl state to check for validity
/// \return True if the pixel corresponding to the state is free, False otherwise
bool isPointValid(cv::Mat& img, const ompl::base::State * state)
{
    double* values = state->as<
      ompl::base::RealVectorStateSpace::StateType>()->values;
    int r = img.rows, c = img.cols;

    cv::Point p((int)(values[0]*c),(int)(values[1]*r));
    
    //std::this_thread::sleep_for(std::chrono::milliseconds(50));

    if(img.at<uchar>(p.y,p.x)==0) { //Black pixel
        return false;
    }

    return true;
}

bool isPointValidSetup(cv::Mat& img, const ompl::base::State * state)
{
    double* values = state->as<
      ompl::base::RealVectorStateSpace::StateType>()->values;
    int r = img.rows, c = img.cols;

    cv::Point p((int)(values[0]*c),(int)(values[1]*r));
    
    if(img.at<uchar>(p.y,p.x)==0) { //Black pixel
        return false;
    }

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


void initializeBeliefModels(cv::Mat& img,
                           const ompl::base::StateSpacePtr space,
                           std::shared_ptr<batching_pomp::cspacebelief::Model<batching_pomp::cspacebelief::BeliefPoint>>& _beliefModel,
                           std::shared_ptr<batching_pomp::cspacebelief::Model<batching_pomp::cspacebelief::BeliefPoint>>& _beliefModel2)
{

  std::random_device rd; 
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis(0.0,1.0);;

  for(unsigned int i = 1u; i <= 1000u; i++) {
    //std::cout<<i<<" : ";
    double x{dis(gen)};
    double y{dis(gen)};

    ompl::base::ScopedState<ompl::base::RealVectorStateSpace> tempState
      = make_state(space,x,y);
    ompl::base::ScopedState<ompl::base::RealVectorStateSpace> tempState2
      = make_state(space,x,y);

    double result;
    if(isPointValidSetup(img,tempState.get())){
      result = 0.0;
    }
    else{
      result = 1.0;
    }

    //std::cout<<"("<<x<<","<<y<<") : "<<result<<std::endl;
    batching_pomp::cspacebelief::BeliefPoint tempBP(tempState.get(),2,result);
    _beliefModel->addPoint(tempBP);
    batching_pomp::cspacebelief::BeliefPoint tempBP2(tempState2.get(),2,result);
    _beliefModel2->addPoint(tempBP2);
  }

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

void displayRoadmap(const ompl::base::StateSpacePtr space,
                    const batching_pomp::Graph& g,
                    const std::string& win_name,
                    cv::Mat& img)
{
  batching_pomp::VertexIter vi, vi_end;
  batching_pomp::EdgeIter ei, ei_end;

  cv::Mat copy = img.clone();
  int r=img.rows,c=img.cols;

  // Display vertices
  for(boost::tie(vi,vi_end)=boost::vertices(g); vi!=vi_end; ++vi)
  {
    ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
      vertex_state(space, g[*vi].v_state->state);
    std::vector<double> coords = vertex_state.reals();
    cv::Point pv((int)(coords[0]*c),(int)(coords[1]*r));
    circle(copy,pv,6,cv::Scalar(0,0,255),-1);
  }

  // Display edges
  for(boost::tie(ei,ei_end)=boost::edges(g); ei!=ei_end; ++ei)
  {
    ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
      v1_state(space, g[boost::source(*ei,g)].v_state->state);
    ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
      v2_state(space, g[boost::target(*ei,g)].v_state->state);
    std::vector<double> v1_coords = v1_state.reals();
    std::vector<double> v2_coords = v2_state.reals();
    cv::Point pv1((int)(v1_coords[0]*c),(int)(v1_coords[1]*r));
    cv::Point pv2((int)(v2_coords[0]*c),(int)(v2_coords[1]*r));

    line(copy,pv1,pv2,cv::Scalar(255,0,0),1);
  }

  //cv::imshow(win_name,copy);
  //cv::waitKey(0);
  cv::imwrite(win_name,copy);
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

double haltonRadiusFun(unsigned int n)
{
  
  auto cardDbl = static_cast<double>(n);

  return 2.5*std::pow(1.0 / cardDbl, 1/ 2.0);

}

double beliefDistanceFunction(
  std::function<double(const ompl::base::State*, const ompl::base::State*)>& _distFun,
  const batching_pomp::cspacebelief::BeliefPoint& _bp1, const batching_pomp::cspacebelief::BeliefPoint& _bp2)
{
  Eigen::VectorXd diff{_bp1.stateValues - _bp2.stateValues};
  return diff.norm();
}

int main(int argc, char* argv[])
{

  po::options_description desc("2d Map Planner Options");

  desc.add_options()
      ("roadmapfile,f",po::value< std::string >()->required(), "Path to graph file")
      ("mapfile,m",po::value< std::string >()->required(),"Environment map image")
      ("start_coords,s",po::value<std::vector<float> >()->multitoken(),"start coordinates")
      ("goal_coords,g",po::value<std::vector<float> >()->multitoken(),"goal coordinates")
      ("num_samples,n",po::value< unsigned int >()->required(),"Number of samples for roadmap")
      ("num_perturbations,p",po::value< unsigned int >()->required(), "Number of times to perturb roadmap samples")
      ("perturb_size,r",po::value<double>()->required(),"Increment for alpha")
      ("dalpha,d",po::value<double>()->required(),"Increment for alpha")
  ;

  // Read arguments
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  std::string map_file(vm["mapfile"].as<std::string>());
  std::string graph_file(vm["roadmapfile"].as<std::string>());
  std::vector<float> start(vm["start_coords"].as<std::vector< float> >());
  std::vector<float> goal(vm["goal_coords"].as<std::vector< float> >());
  unsigned int num_samples(vm["num_samples"].as<unsigned int>());
  unsigned int num_perturbations(vm["num_perturbations"].as<unsigned int>());
  double perturb_size(vm["perturb_size"].as<double>());
  double dalpha(vm["dalpha"].as<double>());

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

  std::shared_ptr<batching_pomp::batching::BatchingManager> dbp;
  dbp.reset(new batching_pomp::batching::DummyBatching(space, graph_file)); 

  std::shared_ptr<batching_pomp::batching::BatchingManager> dbp2;
  dbp2.reset(new batching_pomp::batching::DummyBatching(space, graph_file));

  std::function<double(const ompl::base::State*, const ompl::base::State*)>
    spaceDistFun = std::bind(&ompl::base::RealVectorStateSpace::distance,
                  space->as<ompl::base::RealVectorStateSpace>(),std::placeholders::_1,std::placeholders::_2);
  std::function<double(const batching_pomp::cspacebelief::BeliefPoint& bp1, const batching_pomp::cspacebelief::BeliefPoint& bp2)>
    bpDistFun = std::bind(beliefDistanceFunction,spaceDistFun,std::placeholders::_1,std::placeholders::_2);
  std::shared_ptr<batching_pomp::cspacebelief::Model<batching_pomp::cspacebelief::BeliefPoint>> beliefModel;
  beliefModel.reset(new batching_pomp::cspacebelief::KNNModel(15,0.5,0.2,bpDistFun));

  std::shared_ptr<batching_pomp::cspacebelief::Model<batching_pomp::cspacebelief::BeliefPoint>> beliefModel2;
  beliefModel2.reset(new batching_pomp::cspacebelief::KNNModel(15,0.5,0.2,bpDistFun));

  initializeBeliefModels(img,space,beliefModel,beliefModel2);

  double radius{haltonRadiusFun(num_samples)};
  dbp->setCurrentRadius(radius);

  std::cout<<"Radius of conn. is "<<radius<<std::endl;

  std::shared_ptr< batching_pomp::util::Selector<batching_pomp::Graph>> selectorPtr = 
      std::make_shared<batching_pomp::util::Selector<batching_pomp::Graph>>("normal");

  batching_pomp::BatchingPOMP bpPlanner(si,
                                        dbp,
                                        beliefModel,
                                        selectorPtr,
                                        graph_file,
                                        0.0,
                                        dalpha,
                                        0.05);


  batching_pomp::BatchingPOMP bpPlanner2(si,
                                        dbp2,
                                        beliefModel2,
                                        selectorPtr,
                                        graph_file,
                                        0.0,
                                        dalpha,
                                        0.05);

  bpPlanner.setProblemDefinition(pdef);
  bpPlanner2.setProblemDefinition(pdef);



  batching_pomp::belief_resampling::SGDBasedResampler resampler(bpPlanner,
                                                                num_perturbations,
                                                                perturb_size,
                                                                0.85);
  resampler.setBatchParams(num_samples,radius);
  resampler.updateRoadmap();

  batching_pomp::belief_resampling::SGDBasedResampler resampler2(bpPlanner2,
                                                                0,
                                                                0.0,
                                                                0.85);
  resampler2.setBatchParams(num_samples,radius);
  resampler2.updateRoadmap();

  displayRoadmap(space,*(bpPlanner2.g),"normal.png",img_disp) ;
  displayRoadmap(space,*(bpPlanner.g),"belief.png",img_disp) ;


  ompl::base::PlannerStatus status;
  do
  {
    std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
    status = bpPlanner.solve(ompl::base::plannerNonTerminatingCondition());
    std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
    std::chrono::duration<double> elapsedSeconds{endTime-startTime};
    //planTime += elapsedSeconds.count();

    if(status == ompl::base::PlannerStatus::APPROXIMATE_SOLUTION) {

      auto path = std::dynamic_pointer_cast<ompl::geometric::PathGeometric>(
            pdef->getSolutionPath());

      std::cout<<path->length()<<" "<<bpPlanner.getNumCollChecks()<<std::endl;

      std::size_t pathSize{path -> as< ompl::geometric::PathGeometric >()->getStateCount()};
      state_list.resize(pathSize);

      for(size_t i = 0; i < pathSize; i++) {
        ompl::base::ScopedState<ompl::base::RealVectorStateSpace> st = get_path_state(path -> as< ompl::geometric::PathGeometric >(),i);
        state_list[i] = st.reals();
      }

      //displayPath(img_disp,state_list);
      pdef->clearSolutionPaths();
    }

  }while(status != ompl::base::PlannerStatus::TIMEOUT &&
    status != ompl::base::PlannerStatus::EXACT_SOLUTION);

  std::cout<<"NOW FOR NORMAL ROADMAP"<<std::endl;

  do
  {
    std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
    status = bpPlanner2.solve(ompl::base::plannerNonTerminatingCondition());
    std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
    std::chrono::duration<double> elapsedSeconds{endTime-startTime};
    //planTime += elapsedSeconds.count();

    if(status == ompl::base::PlannerStatus::APPROXIMATE_SOLUTION) {

      auto path = std::dynamic_pointer_cast<ompl::geometric::PathGeometric>(
            pdef->getSolutionPath());

      std::cout<<path->length()<<" "<<bpPlanner2.getNumCollChecks()<<std::endl;

      std::size_t pathSize{path -> as< ompl::geometric::PathGeometric >()->getStateCount()};
      state_list.resize(pathSize);

      for(size_t i = 0; i < pathSize; i++) {
        ompl::base::ScopedState<ompl::base::RealVectorStateSpace> st = get_path_state(path -> as< ompl::geometric::PathGeometric >(),i);
        state_list[i] = st.reals();
      }

      //displayPath(img_disp,state_list);
      pdef->clearSolutionPaths();
    }

  }while(status != ompl::base::PlannerStatus::TIMEOUT &&
    status != ompl::base::PlannerStatus::EXACT_SOLUTION);

  return 0;
}
