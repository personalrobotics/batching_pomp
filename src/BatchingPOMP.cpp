/***********************************************************************
Copyright (c) 2017, Carnegie Mellon University
All rights reserved.
Authors: Shushman Choudhury <shushmanchoudhury@gmail.com>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.
  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#include <functional>
#include <exception>
#include <chrono>
#include <cmath>

#include <ompl/base/goals/GoalState.h>
#include <ompl/base/goals/GoalStates.h>
#include <ompl/geometric/PathGeometric.h>
#include <ompl/util/GeometricEquations.h>
#include <ompl/util/Console.h>
#include <ompl/util/Exception.h>

#include <Eigen/Dense>

#include "batching_pomp/BatchingPOMP.hpp"

namespace batching_pomp{

using batching_pomp::cspacebelief::BeliefPoint;

class ExpWeightMap
{
public:
  typedef boost::readable_property_map_tag category;
  typedef BatchingPOMP::Edge key_type;
  typedef double value_type;
  typedef double reference;
  BatchingPOMP& mPlanner;
  ExpWeightMap(BatchingPOMP& _planner)
  : mPlanner(_planner) {}
};

const double get(const ExpWeightMap& _ewMap, const BatchingPOMP::Edge& e)
{
  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::BLOCKED) {
    return std::numeric_limits<double>::max();
  }

  double alpha{_ewMap.mPlanner.getCurrentAlpha()};

  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::FREE) {
    return (alpha*_ewMap.mPlanner.g[e].distance);
  }

  
  double w_m{-std::log(_ewMap.mPlanner.g[e].probFree)};
  double w_l{_ewMap.mPlanner.g[e].distance};

  double exp_cost{alpha*w_l + (1-alpha)*w_m};
  return exp_cost;
}

class throw_visitor_exception : public std::exception {};
class throw_visitor
{

public:
  BatchingPOMP& mPlanner;
  BatchingPOMP::Vertex mVThrow;
  throw_visitor(BatchingPOMP& _planner,
                BatchingPOMP::Vertex _vThrow)
  : mPlanner(_planner)
  , mVThrow{_vThrow} 
  {}
  inline void initialize_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void discover_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void examine_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g)
  {
    if(u == mVThrow) {
      throw throw_visitor_exception();
    }
  }
  inline void examine_edge(BatchingPOMP::Edge e, const BatchingPOMP::Graph& g)
  {
    mPlanner.computeAndSetEdgeFreeProbability(e);
  }
  inline void edge_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void edge_not_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void black_target(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void finish_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph & g) {}
};


/// For implementing implicit r-neighbour graphs with boost
class neighbours_visitor
{
public:
  BatchingPOMP& mPlanner;
  ompl::NearestNeighbors<BatchingPOMP::Vertex>& mVertexNN;
  double mCurrRadius;
  std::function<double(const BatchingPOMP::Vertex&, const BatchingPOMP::Vertex&)> mVertexDistFun;
  BatchingPOMP::Vertex mVThrow;

  neighbours_visitor(BatchingPOMP& _planner,
                     ompl::NearestNeighbors<BatchingPOMP::Vertex>& _vertexNN, 
                     double _currRadius,
                     BatchingPOMP::Vertex _vThrow)
  : mPlanner(_planner)
  , mVertexNN(_vertexNN)
  , mCurrRadius{_currRadius}
  , mVertexDistFun{mVertexNN.getDistanceFunction()}
  , mVThrow{_vThrow}
  {}
  inline void initialize_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void discover_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void examine_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g)
  {
    if(u == mVThrow) {
      throw throw_visitor_exception();
    }

    std::vector<BatchingPOMP::Vertex> vertexNbrs;
    mVertexNN.nearestR(u,mCurrRadius,vertexNbrs);


    // Now iterate through neighbors and check if edge exists between
    // u and neighbour. If it does not exist, create it AND set its
    // distance (one-time computation) based on provided distance function
    for(BatchingPOMP::Vertex nbr : vertexNbrs){
      // TODO : Find an efficient way to reject itself as a neighbour
      if(nbr == u){
        continue;
      }
      if(!edge(u,nbr,mPlanner.g).second){
        std::pair<BatchingPOMP::Edge,bool> new_edge = add_edge(u,nbr,mPlanner.g);
        mPlanner.g[new_edge.first].distance = mVertexDistFun(source(new_edge.first,mPlanner.g), target(new_edge.first,mPlanner.g));
        mPlanner.g[new_edge.first].blockedStatus = BatchingPOMP::UNKNOWN;
        mPlanner.initializeEdgePoints(new_edge.first);
        mPlanner.computeAndSetEdgeFreeProbability(new_edge.first);
      }
      else{
        // Existing edge - update probability of collision if alpha < 1.0 (measure matters)
        if(mPlanner.getCurrentAlpha() < 1.0) {
          mPlanner.computeAndSetEdgeFreeProbability(edge(u,nbr,mPlanner.g).first);
        }
      }
    }
  }
  inline void examine_edge(BatchingPOMP::Edge e, const BatchingPOMP::Graph& g) {}
  inline void edge_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void edge_not_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void black_target(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void finish_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph & g) {}
};


/// Euclidean distance heuristic for A-star search
template<class Graph, class CostType>
class exp_distance_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
  exp_distance_heuristic(BatchingPOMP& _planner)
  : mPlanner(_planner)
  {}

  CostType operator()(batching_pomp::BatchingPOMP::Vertex u)
  {
    return mPlanner.getCurrentAlpha() * mPlanner.vertexDistFun(u, mPlanner.getGoalVertex());
  }

private:
  BatchingPOMP& mPlanner;
};

template<class Graph, class CostType>
class zero_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
  zero_heuristic(){}

  CostType operator()(batching_pomp::BatchingPOMP::Vertex u)
  {
    return 0.0;
  }
};

/// Computes (symmetric) distance metric between two belief-point instances
/// \param[in] _distFun The distance function between underlying ompl state instances
/// \param[in] _bp1, _bp2 The two belief point instances
/// \return The distance between two belief points
double beliefDistanceFunction(
  std::function<double(const ompl::base::State*, const ompl::base::State*)>& _distFun,
  const BeliefPoint& _bp1, const BeliefPoint& _bp2)
{
  Eigen::VectorXd diff{_bp1.stateValues - _bp2.stateValues};
  return diff.norm();
}

/// Default constructor for CPP-level usage
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
                           std::shared_ptr< batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > _batchingPtr,
                           std::shared_ptr< cspacebelief::Model<BeliefPoint> > _beliefModel,
                           std::unique_ptr< util::Selector<Graph> > _selector,
                           double _increment,
                           double _startGoalRadius,
                           double _pruneThreshold,
                           const std::string& _roadmapFileName)
: ompl::base::Planner(si,"BatchingPOMP")
, mSpace(si->getStateSpace())
, mBatchingPtr{_batchingPtr}
, mBeliefModel{_beliefModel}
, mSelector{std::move(_selector)}
, mCurrentAlpha{0.0}
, mIncrement{_increment}
, mStartGoalRadius{_startGoalRadius}
, mPruneThreshold{_pruneThreshold}
, mIsInitSearchBatch{true}
, mIsPathFound{false}
, mCheckRadius{0.5*mSpace->getLongestValidSegmentLength()}
, mBestPathCost{std::numeric_limits<double>::max()}
, mGraphType{""}
, mBatchingType{""}
, mSelectorType{""}
, mRoadmapName{_roadmapFileName}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mNumLookups{0u}
, mLookupTime{0.0}
, mCollCheckTime{0.0}
, mSearchTime{0.0}
{
  /// Create vertex nearest neighbour manager
  mVertexNN.reset(new ompl::NearestNeighborsGNAT<Vertex>());

  ompl::NearestNeighbors<Vertex>::DistanceFunction distfun(
    [this](const Vertex& a, const Vertex& b)
    {
      return vertexDistFun(a,b);
    });

  mVertexNN->setDistanceFunction(distfun);
}

/// Pure OMPL constructor using parameters
/// Use defaults for objects that cannot be sent with OMPL parameters
/// For others, use zero equivalents
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si)
: ompl::base::Planner(si,"BatchingPOMP")
, mSpace(si->getStateSpace())
, mCurrentAlpha{0.0}
, mIncrement{0.0}
, mStartGoalRadius{0.0}
, mPruneThreshold{DEFAULTPRUNETHRESHOLD}
, mIsInitSearchBatch{true}
, mIsPathFound{false}
, mCheckRadius{0.5*mSpace->getLongestValidSegmentLength()}
, mBestPathCost{std::numeric_limits<double>::max()}
, mGraphType{""}
, mBatchingType{""}
, mSelectorType{""}
, mRoadmapName{""}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mNumLookups{0u}
, mLookupTime{0.0}
, mCollCheckTime{0.0}
, mSearchTime{0.0}
{
  /// Bind distance functions for model manager
  std::function<double(const ompl::base::State*, const ompl::base::State*)>
    spaceDistFun = std::bind(&ompl::base::RealVectorStateSpace::distance,
                  mSpace->as<ompl::base::RealVectorStateSpace>(),std::placeholders::_1,std::placeholders::_2);
  std::function<double(const BeliefPoint& bp1, const BeliefPoint& bp2)>
    bpDistFun = std::bind(beliefDistanceFunction,spaceDistFun,std::placeholders::_1,std::placeholders::_2);

  /// Create model manager with default parameters
  mBeliefModel.reset(new cspacebelief::KNNModel(15,0.5,0.25,bpDistFun));

  /// Create vertex nearest neighbour manager
  mVertexNN.reset(new ompl::NearestNeighborsGNAT<Vertex>());
  ompl::NearestNeighbors<Vertex>::DistanceFunction distfun(
    [this](const Vertex& a, const Vertex& b)
    {
      return vertexDistFun(a,b);
    });

  mVertexNN->setDistanceFunction(distfun);

  /// Define OMPL parameters
  Planner::declareParam<double>("increment", this, &BatchingPOMP::setIncrement, &BatchingPOMP::getIncrement);
  Planner::declareParam<double>("start_goal_radius", this, &BatchingPOMP::setStartGoalRadius, &BatchingPOMP::getStartGoalRadius);
  Planner::declareParam<double>("prune_threshold", this, &BatchingPOMP::setPruneThreshold, &BatchingPOMP::getPruneThreshold);
  Planner::declareParam<std::string>("graph_type", this, &BatchingPOMP::setGraphType, &BatchingPOMP::getGraphType);
  Planner::declareParam<std::string>("batching_type", this, &BatchingPOMP::setBatchingType, &BatchingPOMP::getBatchingType);
  Planner::declareParam<std::string>("selector_type", this, &BatchingPOMP::setSelectorType, &BatchingPOMP::getSelectorType);
  Planner::declareParam<std::string>("roadmap_filename", this, &BatchingPOMP::setRoadmapFileName, &BatchingPOMP::getRoadmapFileName);
}


BatchingPOMP::~BatchingPOMP()
{
}

////////////////////////////////////////////////////////////////////
/// Setters and Getters
double BatchingPOMP::getCurrentAlpha() const
{
  return mCurrentAlpha;
}

double BatchingPOMP::getCurrentBestCost() const
{
  return mBestPathCost;
}

BatchingPOMP::Vertex BatchingPOMP::getStartVertex() const
{
  return mStartVertex;
}

BatchingPOMP::Vertex BatchingPOMP::getGoalVertex() const
{
  return mGoalVertex;
}

bool BatchingPOMP::isInitSearchBatch() const
{
  return mIsInitSearchBatch;
}

double BatchingPOMP::getIncrement() const
{
  return mIncrement;
}

void BatchingPOMP::setIncrement(double _increment)
{
  OMPL_INFORM("Increment set!");
  mIncrement = _increment;
}

double BatchingPOMP::getStartGoalRadius() const
{
  return mStartGoalRadius;
}

void BatchingPOMP::setStartGoalRadius(double _startGoalRadius)
{
  mStartGoalRadius = _startGoalRadius;
}

double BatchingPOMP::getPruneThreshold() const
{
  return mPruneThreshold;
}

void BatchingPOMP::setPruneThreshold(double _pruneThreshold)
{
  mPruneThreshold = _pruneThreshold;
}


std::string BatchingPOMP::getGraphType() const
{
  return mGraphType;
}

void BatchingPOMP::setGraphType(const std::string& _graphType)
{
  mGraphType = _graphType;
}

std::string BatchingPOMP::getBatchingType() const
{
  return mBatchingType;
}

void BatchingPOMP::setBatchingType(const std::string& _batchingType)
{
  mBatchingType = _batchingType;
}

std::string BatchingPOMP::getSelectorType() const
{
  return mSelectorType;
}

void BatchingPOMP::setSelectorType(const std::string& _selectorType)
{
  mSelectorType = _selectorType;
}

std::string BatchingPOMP::getRoadmapFileName() const
{
  return mRoadmapName;
}

void BatchingPOMP::setRoadmapFileName(const std::string& _roadmapFileName)
{
  mRoadmapName = _roadmapFileName;
}


/// Public helper methods
double BatchingPOMP::vertexDistFun(const Vertex& u, const Vertex& v) const
{
  return mSpace->distance(g[u].v_state->state, g[v].v_state->state);
}


void BatchingPOMP::initializeEdgePoints(const Edge& e)
{
  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;

  unsigned int nStates = static_cast<unsigned int>(std::floor(g[e].distance / (2.0*mCheckRadius)));
  
  if(nStates < 2u) {
    nStates = 2u; // Just start and goal
  }

  g[e].edgeStates.resize(nStates);

  for(unsigned int i = 0; i < nStates; i++)
  {
    g[e].edgeStates[i].reset(new StateCon(mSpace));
  }

  const std::vector< std::pair<int,int> > & order = mBisectPermObj.get(nStates);

  for(unsigned int i = 0; i < nStates; i++)
  {
    mSpace->interpolate(startState, endState,
      1.0*(1+order[i].first)/(nStates+1), g[e].edgeStates[i]->state);
  }
}


double BatchingPOMP::computeAndSetEdgeFreeProbability(const Edge& e)
{
  /// March along edge states with some jump factor
  /// And estimate the probability of collision of each

  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;

  auto nStates = g[e].edgeStates.size();
  unsigned int stepSize{sqrt(nStates)};
  double result{1.0};

  for(unsigned int i = 0; i < nStates; i+=stepSize)
  {
    BeliefPoint query(g[e].edgeStates[i]->state, mSpace->getDimension(), -1.0);

    /// Timing and counting block for model estimation
    mNumLookups++;
    std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
    double collProb{mBeliefModel->estimate(query)};
    std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
    std::chrono::duration<double> elapsedSeconds{endTime-startTime};
    mLookupTime += elapsedSeconds.count();

    result *= (1.0 - collProb);
  }

  g[e].probFree = result;

  return result;
}


bool BatchingPOMP::checkAndSetEdgeBlocked(const BatchingPOMP::Edge& e)
{
  /// March along edge states with highest resolution
  mNumEdgeChecks++;

  auto validityChecker = si_->getStateValidityChecker();
  
  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;

  auto nStates = g[e].edgeStates.size();

  for(unsigned int i = 1; i < nStates-1; i++)
  {
    ompl::base::ScopedState<> sstate(mSpace);
    sstate = g[e].edgeStates[i]->state;

    /// Counting and timing block for collision checks
    mNumCollChecks++;
    std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
    bool checkResult{validityChecker->isValid(g[e].edgeStates[i]->state)};
    std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
    std::chrono::duration<double> elapsedSeconds{endTime-startTime};
    mCollCheckTime += elapsedSeconds.count();
    
    if(checkResult == false) {

      BeliefPoint toAdd(g[e].edgeStates[i]->state,mSpace->getDimension(),1.0);
      mBeliefModel->addPoint(toAdd);

      /// Update edge properties to reflect blocked
      g[e].blockedStatus = BatchingPOMP::BLOCKED;
      g[e].distance = std::numeric_limits<double>::max();
      g[e].probFree = 0.0;
      return true;
    }
    else {
      BeliefPoint toAdd(g[e].edgeStates[i]->state,mSpace->getDimension(),0.0);
      mBeliefModel->addPoint(toAdd);
    }
  }

  /// Update edge properties to reflect free
  g[e].blockedStatus = BatchingPOMP::FREE;
  g[e].probFree = 1.0;
  return false;
}


/// Private helper methods
double BatchingPOMP::haltonRadiusFun(unsigned int n) const
{
  
  auto dimDbl = static_cast<double>(mSpace->getDimension());
  auto cardDbl = static_cast<double>(n);

  return 2.5*std::pow(1.0 / cardDbl, 1/ dimDbl);

}

double BatchingPOMP::rggRadiusFun(unsigned int n) const
{
  /// Lifted from BIT* - need to find a way to use directly
  auto dimDbl = static_cast<double>(mSpace->getDimension());

  double approximationMeasure{si_->getSpaceMeasure()};
  if (!std::isfinite(approximationMeasure)) {
    throw ompl::Exception("Measure of space is unbounded!");
  }

  double minRggR{
    2.0 * std::pow((1.0 + 1.0 / dimDbl) * (approximationMeasure / ompl::unitNBallMeasure(si_->getStateDimension())), 1.0/dimDbl)};

  auto cardDbl = static_cast<double>(n);

  return minRggR * std::pow(std::log(cardDbl) / cardDbl, 1 / dimDbl);
}


bool BatchingPOMP::checkAndUpdatePathBlocked(const std::vector<Edge>& _ePath)
{
  std::vector<Edge> selectedEPath = mSelector->selectEdges(g,_ePath);

  for(Edge e : selectedEPath)
  {
    if(g[e].blockedStatus == BatchingPOMP::UNKNOWN) {
      
      if(checkAndSetEdgeBlocked(e)) {
        return true;
      }
    }
  }
  return false;
}


bool BatchingPOMP::isVertexInadmissible(const ompl::base::State* vState) const
{
  if(mBestPathCost >= std::numeric_limits<double>::max()) {
    return false;
  }

  double bestCostThroughVertex
    {mSpace->distance(g[mStartVertex].v_state->state,vState)
      + mSpace->distance(g[mGoalVertex].v_state->state,vState)};

  return (bestCostThroughVertex > mBestPathCost);
}


bool BatchingPOMP::vertexPruneFunction(const ompl::base::State* vState) const
{
  auto validityChecker = si_->getStateValidityChecker();

  return (isVertexInadmissible(vState)==true || validityChecker->isValid(vState)==false);

}

double BatchingPOMP::getPathDistance(const std::vector<Edge>& _ePath) const
{
  double pathDistance{0.0};
  for(auto e : _ePath)
  {
    pathDistance += g[e].distance;
  }
  return pathDistance;
}

void BatchingPOMP::setup()
{

  // TODO - Make scaling params

  Planner::setup();

  /// Radius function type for Edge or Hybrid Batching
  if(mGraphType == "halton") {
    mRadiusFun = std::bind(&BatchingPOMP::haltonRadiusFun,this,std::placeholders::_1);
  }
  else if (mGraphType == "rgg") {
    mRadiusFun = std::bind(&BatchingPOMP::rggRadiusFun,this,std::placeholders::_1);
  }
  else{
    if(mBatchingType=="edge" || mBatchingType=="hybrid") {
      throw ompl::Exception("Invalid graph type specified! Graph type needed for hybrid or edge batching!");
    }
  }

  if(mRoadmapName == "") {
    throw ompl::Exception("Roadmap name must be set before creating batching manager!");
  }

  /// Now create batching pointer with default parameters
  if(mBatchingType == "vertex") {
    unsigned int initNumVertices{100u};
    double vertInflFactor{2.0};
    
    std::shared_ptr<batching::VertexBatching<Graph,VPStateMap,StateCon,EPDistanceMap>> vbp 
       = std::make_shared<batching::VertexBatching<Graph,VPStateMap,StateCon,EPDistanceMap>>
        (mSpace, get(&VProps::v_state,full_g), mRoadmapName, full_g, g, initNumVertices, vertInflFactor);

    mBatchingPtr = std::static_pointer_cast<batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>>
                  (std::move(vbp));
  }
  else if(mBatchingType == "edge") {
    auto dimDbl = static_cast<double>(mSpace->getDimension());
    double radiusInflFactor{std::pow(2.0,1/dimDbl)};
    double maxRadius{mSpace->getMaximumExtent()};

    std::shared_ptr<batching::EdgeBatching<Graph,VPStateMap,StateCon,EPDistanceMap>> ebp
      = std::make_shared<batching::EdgeBatching<Graph,VPStateMap,StateCon,EPDistanceMap>>
        (mSpace, get(&VProps::v_state,full_g), mRoadmapName, full_g, g, radiusInflFactor, mRadiusFun, maxRadius);
    
    mBatchingPtr = std::static_pointer_cast<batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>>
                  (std::move(ebp));
  }
  else if(mBatchingType == "hybrid") {
    unsigned int initNumVertices{100u};
    double vertInflFactor{2.0};
    auto dimDbl = static_cast<double>(mSpace->getDimension());
    double radiusInflFactor{std::pow(2.0,1/dimDbl)};
    double maxRadius{mSpace->getMaximumExtent()};

    std::shared_ptr<batching::HybridBatching<Graph,VPStateMap,StateCon,EPDistanceMap>> hbp
     = std::make_shared<batching::HybridBatching<Graph,VPStateMap,StateCon,EPDistanceMap>>
       (mSpace, get(&VProps::v_state,full_g), mRoadmapName, full_g, g, initNumVertices, vertInflFactor, radiusInflFactor, mRadiusFun, maxRadius);

    mBatchingPtr = std::static_pointer_cast<batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>>
                  (std::move(hbp));

  }
  else if(mBatchingType == "single") {

    std::shared_ptr<batching::SingleBatching<Graph,VPStateMap,StateCon,EPDistanceMap>> sbp 
     = std::make_shared<batching::SingleBatching<Graph,VPStateMap,StateCon,EPDistanceMap>>
        (mSpace, get(&VProps::v_state,full_g), get(&EProps::distance,full_g), mRoadmapName, full_g, g);
    
    mBatchingPtr = std::static_pointer_cast<batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>>
                  (std::move(sbp));

    // Setup all edges
    EdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=edges(g); ei!=ei_end; ++ei)
    {
      initializeEdgePoints(*ei);
    }
  }
  else {
    throw ompl::Exception("Invalid batching type specified - "+mBatchingType+"!");
  }

  /// Create selector with the type specified
  if(mSelectorType == "") {
    throw ompl::Exception("Selector type must be set before creation");
  }
  mSelector.reset(new batching_pomp::util::Selector<Graph>(mSelectorType));
}


void BatchingPOMP::setProblemDefinition(
  const ompl::base::ProblemDefinitionPtr & pdef)
{
  ompl::base::Planner::setProblemDefinition(pdef);

  StateConPtr startState{std::make_shared<StateCon>(mSpace)};
  mSpace->copyState(startState->state, pdef->getStartState(0));

  ompl::base::GoalPtr goal{pdef->getGoal()};
  ompl::base::GoalState* goalStatePtr{goal->as<ompl::base::GoalState>()};
  StateConPtr goalState{std::make_shared<StateCon>(mSpace)};
  mSpace->copyState(goalState->state, goalStatePtr->getState());

  auto validityChecker = si_->getStateValidityChecker();

  /// Just add start and goal to empty roadmap and to vertex ptr
  /// Batching will take care of adding neighbours and so on
  if(!validityChecker->isValid(startState->state)) {
    throw ompl::Exception("Start configuration is in collision!");
  }
  mStartVertex = boost::add_vertex(g);
  g[mStartVertex].v_state = startState;
  mVertexNN->add(mStartVertex);

  if(!validityChecker->isValid(goalState->state)) {
    throw ompl::Exception("Goal configuration is in collision!");
  }
  mGoalVertex = boost::add_vertex(g);
  g[mGoalVertex].v_state = goalState;
  mVertexNN->add(mGoalVertex);

  /// For single batch, add start and goal to roadmap with startGoalRadius
  if(mBatchingType == "single") {

    VertexIter vi, vi_end;
    for(boost::tie(vi,vi_end)=vertices(g); vi!=vi_end; ++vi) {
      if(*vi == mStartVertex || *vi== mGoalVertex) {
        continue;
      }
      double startDist = vertexDistFun(mStartVertex,*vi);
      if(startDist < mStartGoalRadius) {
        std::pair<Edge, bool> new_edge = add_edge(*vi,mStartVertex,g);
        g[new_edge.first].blockedStatus = UNKNOWN;
        g[new_edge.first].distance = startDist;
        initializeEdgePoints(new_edge.first);
      }

      double goalDist = vertexDistFun(mGoalVertex,*vi);
      if(goalDist < mStartGoalRadius) {
        std::pair<Edge, bool> new_edge = add_edge(*vi,mGoalVertex,g);
        g[new_edge.first].blockedStatus = UNKNOWN;
        g[new_edge.first].distance = goalDist;
        initializeEdgePoints(new_edge.first);
      }
    }
  }
}


ompl::base::PlannerStatus BatchingPOMP::solve(
  const ompl::base::PlannerTerminationCondition & ptc)
{
  
  double currSolnCost{std::numeric_limits<double>::max()};
  mIsPathFound = false;
  std::vector<Edge> ePath;

  while(currSolnCost >= mBestPathCost)
  {

    if(mBatchingPtr->isExhausted() && mIsInitSearchBatch)
    {
      if(mBestPathCost < std::numeric_limits<double>::max()) {
        /// Non trivial best solution found
        OMPL_INFORM("All batches exhausted - current solution is the best!");
        return ompl::base::PlannerStatus::EXACT_SOLUTION;
      }
      else {
        /// No solution found
        OMPL_INFORM("All batches exhausted - no solution found!");
        return ompl::base::PlannerStatus::TIMEOUT;
      }
    }

    if(ptc == true){
      OMPL_INFORM("Planner termination condition satisfied!");
      return ompl::base::PlannerStatus::TIMEOUT;
    }

    /// If it is the start of a new batch, get the next Batch
    /// And reset current inflation factor
    if(mIsInitSearchBatch) {
      mCurrentAlpha = 0.0;
      
      std::function<bool(const ompl::base::State*)> pruneFunction = 
        std::bind(&BatchingPOMP::vertexPruneFunction, this,std::placeholders::_1);
      mBatchingPtr->nextBatch(pruneFunction, *mVertexNN);
    }

    /// Now do the search
    std::map<Vertex,Vertex> startPreds;
    std::map<Vertex,double> startDist;
    std::map<Vertex,double> startFValue;
    std::map<Vertex,boost::default_color_type> colorMap;

    try
    {
      mNumSearches++;
      std::unique_ptr<boost::astar_heuristic<Graph, double>> heuristicFn;

      if(mIsInitSearchBatch) {
        heuristicFn.reset(new zero_heuristic<Graph,double>());
      }
      else {
        heuristicFn.reset(new exp_distance_heuristic<Graph,double>(*this));
      }

      mIsInitSearchBatch = false; // Either way, make it false

      /// TODO : Assign to a generic visitor object?
      std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
      if(mBatchingType == "single") {
        boost::astar_search(
          g,
          mStartVertex,
          *heuristicFn,
          throw_visitor(*this, mGoalVertex),
          boost::make_assoc_property_map(startPreds),
          boost::make_assoc_property_map(startFValue),
          boost::make_assoc_property_map(startDist),
          ExpWeightMap(*this),
          boost::get(boost::vertex_index, g),
          boost::make_assoc_property_map(colorMap),
          std::less<double>(),
          boost::closed_plus<double>(std::numeric_limits<double>::max()),
          std::numeric_limits<double>::max(),
          double()
        );
      }
      else {
        boost::astar_search(
          g,
          mStartVertex,
          *heuristicFn,
          neighbours_visitor(*this, *mVertexNN, mBatchingPtr->getCurrentRadius(), mGoalVertex),
          boost::make_assoc_property_map(startPreds),
          boost::make_assoc_property_map(startFValue),
          boost::make_assoc_property_map(startDist),
          ExpWeightMap(*this),
          boost::get(boost::vertex_index, g),
          boost::make_assoc_property_map(colorMap),
          std::less<double>(),
          boost::closed_plus<double>(std::numeric_limits<double>::max()),
          std::numeric_limits<double>::max(),
          double()
        );
      }
      std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
      std::chrono::duration<double> elapsedSeconds{endTime-startTime};
      mSearchTime += elapsedSeconds.count();
    }
    catch (const throw_visitor_exception & ex)
    {
    }
    
    if(startDist[mGoalVertex] == std::numeric_limits<double>::max()) {
      /// Did not find a path this time
      OMPL_INFORM("No new feasible path found in this batch!");
      if(mBatchingPtr->isExhausted()) {
        if(mBestPathCost < std::numeric_limits<double>::max()) {
          OMPL_INFORM("No further solutions on roadmap");
          return ompl::base::PlannerStatus::EXACT_SOLUTION;
        }
        else {
          OMPL_INFORM("No solution found on roadmap");
          return ompl::base::PlannerStatus::TIMEOUT;
        }
      }
      else {
        // Trigger a new batch
        mIsInitSearchBatch = true;
        if(mBestPathCost < std::numeric_limits<double>::max()) {
          return ompl::base::PlannerStatus::APPROXIMATE_SOLUTION;
        }
        else {
          // ABORT so solve() can be run again
          return ompl::base::PlannerStatus::ABORT;
        }
      }

    }

    /// Retrieve path and evaluate
    Vertex vWalk{mGoalVertex};

    ePath.clear();
    while(vWalk != mStartVertex)
    {
      Vertex vPred{startPreds[vWalk]};
      std::pair<Edge, bool> edgePair{boost::edge(vPred, vWalk, g)};
      if(!edgePair.second) {
        throw ompl::Exception("Error! Edge present during search no longer exists in graph!");
      }
      ePath.push_back(edgePair.first);
      vWalk = vPred;
    }
    if(ePath.size() > 1) {
      std::reverse(ePath.begin(), ePath.end());
    }
    /// If path is free, set current cost to it, increment alpha
    bool pathBlocked{checkAndUpdatePathBlocked(ePath)};

    if(!pathBlocked) {
      mIsPathFound = true;
      currSolnCost = getPathDistance(ePath);

      if(mCurrentAlpha >= 1.0) {
        mIsInitSearchBatch = true;
      }
      else {
        mCurrentAlpha = std::min(mCurrentAlpha+mIncrement, 1.0);
      }
    }
  }


  /// Has found an improved feasible solution -> now report
  double improvementFactor{mBestPathCost/currSolnCost};
  mBestPathCost = currSolnCost;
  mCurrBestPath = ePath;

  mBatchingPtr->updateWithNewSolutionCost(mBestPathCost);

  ompl::geometric::PathGeometric* path = new ompl::geometric::PathGeometric(si_);
  VertexIndexMap vertex_id = get(boost::vertex_index,g);

  path->append(g[mStartVertex].v_state->state);

  for(auto e : mCurrBestPath)
  {
    path->append(g[target(e,g)].v_state->state);
  }

  pdef_->addSolutionPath(ompl::base::PathPtr(path));

  /// Prune vertices using current solution cost if sufficiently improved
  if(improvementFactor > mPruneThreshold) {
    std::function<bool(const ompl::base::State* vState)> pruneFunction = 
        std::bind(&BatchingPOMP::isVertexInadmissible, this,std::placeholders::_1);
    mBatchingPtr->pruneVertices(pruneFunction,*mVertexNN);
  }

  return ompl::base::PlannerStatus::APPROXIMATE_SOLUTION;

}


} // namespace batching_pomp