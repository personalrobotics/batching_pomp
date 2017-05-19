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

#include "batching_pomp/include/BatchingPOMP.hpp"

namespace batching_pomp{

using batching_pomp::cspacebelief::BeliefPoint;

/// Default constructor for CPP-level usage
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
                           std::unique_ptr< BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > _batchingPtr,
                           std::unique_ptr< Model<BeliefPoint> > _beliefModel,
                           std::unique_ptr< Selector<Graph> > _selector,
                           double _increment,
                           double _startGoalRadius,
                           const std::string& _roadmapFileName)
: ompl::base::Planner(si,"BatchingPOMP")
, mSpace(si->getStateSpace()),
, mBatchingPtr{std::move(_batchingPtr)}
, mBeliefModel{std::move(_beliefModel)}
, mSelector{std::move(_selector)}
, mCurrentAlpha{0.0}
, mIncrement{_increment}
, mStartGoalRadius{_startGoalRadius}
, mIsInitSearchBatch{true}
, mIsPathfound{false}
, mCheckRadius{0.5*space->getLongestValidSegmentLength()}
, mBestPathCost{std::numeric_limits<double>::max()}
, mGraphType{""}
, mBatchingType{""}
, mSelectorType{""}
, mRoadmapName{_roadmapFileName}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mLookupTime{0.0}
{
  /// Create vertex nearest neighbour manager
  mVertexNN.reset(new ompl::NearestNeighborsGNAT<Vertex>());
  mVertexNN->setDistanceFunction(vertexDistFun);
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
, mIsInitSearchBatch{true}
, mIsPathfound{false}
, mCheckRadius{0.5*space->getLongestValidSegmentLength()}
, mBestPathCost{std::numeric_limits<double>::max()}
, mGraphType{""}
, mBatchingType{""}
, mSelectorType{""}
, mRoadmapName{""}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mLookupTime{0.0}
{
  /// Bind distance functions for model manager
  std::function<double(const ompl::base::State*, const ompl::base::State*)>
    spaceDistFun = std::bind(&ompl::base::RealVectorStateSpace::distance,
                  space->as<ompl::base::RealVectorStateSpace>(),std::placeholders::_1,std::placeholders::_2);
  std::function<double(const BeliefPoint& bp1, const BeliefPoint& bp2)>
    bpDistFun = std::bind(batching_pomp::cspacebelief::beliefDistanceFunction,spaceDistFun,std::placeholders::_1,std::placeholders::_2);

  /// Create model manager with default parameters
  mBeliefModel.reset(new batching_pomp::cspacebelief::KNNModel(15,2.0,0.5,0.25,bpDistFun));

  /// Create vertex nearest neighbour manager
  mVertexNN.reset(new ompl::NearestNeighborsGNAT<Vertex>());
  mVertexNN->setDistanceFunction(vertexDistFun);

  /// Define OMPL parameters
  Planner::declareParam<double>("increment", this, &BatchingPOMP::setIncrement, &BatchingPOMP::getIncrement);
  Planner::declareParam<double>("start_goal_radius", this, &BatchingPOMP::setStartGoalRadius, &BatchingPOMP::getStartGoalRadius);
  Planner::declareParam<std::string>("graph_type", this, &BatchingPOMP::setGraphType, &BatchingPOMP::getGraphType);
  Planner::declareParam<std::string>("batching_type", this, &BatchingPOMP::setBatchingType, &BatchingPOMP::getBatchingType);
  Planner::declareParam<std::string>("selector_type", this, &BatchingPOMP::setSelectorType, &BatchingPOMP::getSelectorType);
  Planner::declareParam<std::string>("roadmap_filename", this, &BatchingPOMP::setRoadmapFileName, &BatchingPOMP::getRoadmapFileName);

  /// Create objects that depend on parameters

  /// Radius function type for Edge or Hybrid Batching
  if(mGraphType == "" && (mBatchingType=="edge" || mBatchingType=="hybrid")) {
    throw std::runtime_error("Need a graph type to get radius function!");
  }
  if(mGraphType == "halton") {
    mRadiusFun = haltonRadiusFun;
  }
  else if (mGraphType == "rgg") {
    mRadiusFun = rggRadiusFun;
  }
  else{
    throw std::invalid_argument("Invalid graph type specified - "+mGraphType+"!");
  }

  if(mRoadmapName == "") {
    throw std::runtime_error("Roadmap name must be set before creating batching manager!");
  }

  /// Now create batching pointer with default parameters
  if(mBatchingType == "vertex") {
    unsigned int initNumVertices = 100u;
    double vertInflFactor = 2.0;
    
    mBatchingPtr.reset(
      new batching_pomp::batching::VertexBatching<Graph,VPStateMap,StateCon>
        (mSpace, get(&VProps::v_state,g), mRoadmapName, g, initNumVertices, vertInflFactor) );
  }
  else if(mBatchingType == "edge") {
    auto dimDbl = static_cast<double>(mSpace->getDimension());
    double radiusInflFactor = std::pow(2.0,1/dimDbl);
    double maxRadius = mSpace->getMaximumExtent();
    
    mBatchingPtr.reset(
      new batching_pomp::batching::EdgeBatching<Graph,VPStateMap,StateCon>
        (mSpace, get(&VProps::v_state,g), mRoadmapName, g, radiusInflFactor, mRadiusFun, maxRadius) );
  }
  else if(mBatchingType == "hybrid") {
    unsigned int initNumVertices = 100u;
    double vertInflFactor = 2.0;
    auto dimDbl = static_cast<double>(mSpace->getDimension());
    double radiusInflFactor = std::pow(2.0,1/dimDbl);
    double maxRadius = mSpace->getMaximumExtent();

    mBatchingPtr.reset(
      new batching_pomp::batching::HybridBatching<Graph,VPStateMap,StateCon>
        (mSpace, get(&VProps::v_state,g), mRoadmapName, g, initNumVertices, vertInflFactor, radiusInflFactor, mRadiusFun, maxRadius) );
  }
  else if(mBatchingType == "single") {
    mBatchingPtr.reset(
      new batching_pomp::batching::SingleBatching<Graph,VPStateMap,StateCon,EPDistanceMap>
       (mSpace, get(&VProps::v_state,g), get(&EProps::distance,g), mRoadmapName, g ) );
  }
  else {
    throw std::runtime_error("Invalid batching type specified - "+mBatchingType+"!");
  }

  /// Create selector with the type specified
  if(mSelector == "") {
    throw std::runtime_error("Selector type must be set before creation");
  }
  mSelector.reset(new batching_pomp::util::Selector(mSelectorType));
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

Vertex BatchingPOMP::getStartVertex() const
{
  return mStartVertex;
}

Vertex BatchingPOMP::getGoalVertex() const
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

// double BatchingPOMP::getCheckRadius() const
// {
//   return mCheckRadius;
// }

// void BatchingPOMP::setCheckRadius(double _checkRadius)
// {
//   mCheckRadius = _checkRadius;
// }

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


double BatchingPOMP::computeAndSetEdgeFreeProbability(const Edge& e)
{
  /// March along edge states with some jump factor
  /// And estimate the probability of collision of each

  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;

  unsigned int nStates{std::floor(g[e].distance / (2.0*mCheckRadius))};
  unsigned int stepSize{5};
  double result{1.0};

  for(unsigned int i = 0; i < nStates; i+=stepSize)
  {
    StateConPtr tempState{std::make_shared<StateCon>(mSpace)};
    mSpace->interpolate(startState, endState,
      1.0*(1+i)/(nStates+1), tempState->state);

    BeliefPoint query(tempState->state, -1.0);

    /// Timing block for model estimation
    std::chrono::time_point<std::chrono::system_clock> startTime{std::chrono::system_clock::now()};
    double collProb{mBeliefModel->estimate(query)};
    std::chrono::time_point<std::chrono::system_clock> endTime{std::chrono::system_clock::now()};
    std::chrono::duration<double> elapsedSeconds{endTime-startTime};
    mLookupTime += elapsedSeconds.count();

    result *= (1.0 - coll_prob);
  }

  g[e].probFree = result;
  return result;
}


bool checkAndSetEdgeBlocked(const Edge& e)
{
  /// March along edge states with highest resolution
  
  mNumEdgeChecks++;
  addAffectedEdges(e);

  auto validityChecker = si_->getStateValidityChecker();
  
  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;
  unsigned int nStates{std::floor(g[e].distance / (2.0*mCheckRadius))};

  const std::vector< std::pair<int,int> > & order = bisectPermObj_.get(nStates);

  for(unsigned int i = 0; i < nStates; i++)
  {
    StateConPtr tempState{std::make_shared<StateCon>(mSpace)};
    mSpace->interpolate(startState, endState,
      1.0*(1+order[ui].first)/(nStates+1), tempState->state);

    mNumCollChecks++;

    /// Check and add to belief model
    if(validityChecker->isValid(tempState->state) == false) {
      BeliefPoint toAdd(tempState->state,1.0);
      mBeliefModel->addPoint(toAdd);

      /// Update edge properties to reflect blocked
      g[e].blockedStatus = BatchingPOMP::BLOCKED;
      g[e].distance = std::numeric_limits<double>::max();
      g[e].probFree = std::numeric_limits<double>::max();
      return true;
    }
    else {
      BeliefPoint toAdd(tempState->state,0.0);
      mBeliefModel->addPoint(toAdd);
    }
  }

  /// Update edge properties to reflect free
  g[e].blockedStatus = BatchingPOMP::FREE;
  g[e].probFree = 0.0;
  return false;
}


/// Private helper methods
double BatchingPOMP::haltonRadiusFun(unsigned int n) const
{
  
  auto dimDbl = static_cast<double>(mSpace->getDimension());
  auto cardDbl = static_cast<double>(n);

  return std::pow(1.0 / cardDbl, 1/ dimDbl);

}

double BatchingPOMP::rggRadiusFun(unsigned int n) const
{
  /// Lifted from BIT* - need to find a way to use directly
  auto dimDbl = static_cast<double>(mSpace->getDimension());

  double approximationMeasure{si_->getSpaceMeasure()};
  if (!std::isfinite(approximationMeasure_)) {
    throw ompl::Exception("Measure of space is unbounded!");
  }

  double minRggR{
    2.0 * std::pow((1.0 + 1.0 / dimDbl) * (approximationMeasure_ / unitNBallMeasure(si_->getStateDimension())), 1.0/dimDbl)};

  auto cardDbl = static_cast<double>(n);

  return minRggR * std::pow(std::log(cardDbl) / cardDbl, 1 / dimDbl);
}

void BatchingPOMP::updateAffectedEdgeWeights()
{
  for(auto e : mEdgesToUpdate)
  {
    if(g[e].blockedStatus == BatchingPOMP::UNKNOWN) {
      computeAndSetEdgeFreeProbability(e);
    }
  }

  mEdgesToUpdate.clear();
}


void BatchingPOMP::addAffectedEdges(const Edge& e)
{
  // TODO : Change this to be based on model radius (and what edges are currently in there)
  // I.E you'll want to lookup mVertexNN and add edge if it exists

  /// For now, just add all out_edges and out_edges from neigbours
  Vertex u{source(e,g)};
  Vertex v{target(e,g)};

  OutEdgeIter ei;
  OutEdgeIter ei_end;

  for(boost::tie(ei,ei_end)=out_edges(u,g); ei!=ei_end; ++ei)
  {
    if(g[*ei].blockedStatus == BatchingPOMP::UNKNOWN) {

      /// Add edge if not already there
      if(find(mEdgesToUpdate.begin(), mEdgesToUpdate.end(), *ei) != mEdgesToUpdate.end()) {
        mEdgesToUpdate.push_back(*ei);
      }

      Vertex nbr = target(*ei,g);
      OutEdgeIter ei_in,ei_in_end;

      for(boost::tie(ei_in,ei_in_end)=out_edges(nbr,g); ei_in!=ei_in_end; ++ei)
      {
        Vertex nbr_nbr = target(*ei_in,g);
        if(find(mEdgesToUpdate.begin(), mEdgesToUpdate.end(), *ei_in) != mEdgesToUpdate.end()) {
          mEdgesToUpdate.push_back(*ei_in);
        }
      }

    }
  }

  /// Do the same for target
  for(boost::tie(ei,ei_end)=out_edges(v,g); ei!=ei_end; ++ei)
  {
    if(g[*ei].blockedStatus == BatchingPOMP::UNKNOWN) {
      /// Add edge if not already there
      if(find(mEdgesToUpdate.begin(), mEdgesToUpdate.end(), *ei) != mEdgesToUpdate.end()) {
        mEdgesToUpdate.push_back(*ei);
      }

      Vertex nbr = target(*ei,g);
      OutEdgeIter ei_in,ei_in_end;

      for(boost::tie(ei_in,ei_in_end)=out_edges(nbr,g); ei_in!=ei_in_end; ++ei)
      {
        Vertex nbr_nbr = target(*ei_in,g);
        if(find(mEdgesToUpdate.begin(), mEdgesToUpdate.end(), *ei_in) != mEdgesToUpdate.end()) {
          mEdgesToUpdate.push_back(*ei_in);
        }
      }

    }
  }

}


bool BatchingPOMP::checkAndUpdatePathBlocked(const std::vector<Edge>& _ePath)
{
  std::vector<Edge> selectedEPath = mSelector->selectEdges(g,_ePath);

  for(auto e : selectedEPath)
  {
    if(g[e].blockedStatus == BatchingPOMP::UNKNOWN) {
      
      if(checkAndSetEdgeBlocked(e)) {
        return true;
      }
    }
  }
  return false;
}


bool BatchingPOMP::isVertexAdmissible(const Vertex& v) const
{
  if(mBestPathCost == std::numeric_limits<double>::max()) {
    return true;
  }

  double bestCostThroughVertex{vertexDistFun(mStartVertex,v) + vertexDistFun(v,mGoalVertex)};

  return (bestCostThroughVertex < mBestPathCost);
}


bool BatchingPOMP::vertexPruneFunction(const Vertex& v) const
{
  auto validityChecker = si_->getStateValidityChecker();

  return (isVertexAdmissible(v)==false || validityChecker->isValid(g[v].v_state->state)==false);

}

double getPathDistance(const std::vector<Edge>& _ePath) const
{
  double pathDistance{0.0};
  for(auto e : _ePath)
  {
    pathDistance += g[e].distance;
  }
  return pathDistance;
}


void BatchingPOMP::setProblemDefinition(
  const ompl::base::ProblemDefinitionPtr & pdef)
{
  ompl::base::Planner::setProblemDefinition(pdef);

  StateConPtr startState{std::make_shared<StateCon>(mSpace)};
  mSpace->copyState(startState->state, pdef->getStartState(0));

  ompl::base::GoalPtr goal = pdef->getGoal();
  ompl::base::GoalState* goalStatePtr = goal->as<ompl::base::GoalState>();
  StateConPtr goalState{std::make_shared<StateCon>(mSpace)};
  mSpace->copyState(goalState->state, goalStatePtr->getState());

  auto validityChecker = si_->getStateValidityChecker();

  /// Just add start and goal to empty roadmap and to vertex ptr
  /// Batching will take care of adding neighbours and so on
  if(!validityChecker->isValid(startState->state)) {
    throw std::runtime_error("Start configuration is in collision!");
  }
  mStartVertex = boost::add_vertex(g);
  g[mStartVertex].v_state = startState;
  mVertexNN->add(mStartVertex);

  if(!validityChecker->isValid(goalState->state)) {
    throw std::runtime_error("Start configuration is in collision!");
  }
  mGoalVertex = boost::add_vertex(g);
  g[mGoalVertex].v_state = goalState;
  mVertexNN->add(mGoalVertex);
}


ompl::base::PlannerStatus BatchingPOMP::solve(
  const ompl::base::PlannerTerminationCondition & ptc)
{
  
  double currSolnCost{std::numeric_limits<double>::max()};
  mIsPathfound = false;
  std::vector<Vertex> ePath;

  while(currSolnCost >= mBestPathCost)
  {

    if(mBatchingPtr->isExhausted() && mIsInitSearchBatch)
    {
      if(mBestPathCost < std::numeric_limits<double>::max()) {
        /// Non trivial best solution found
        OMPL_INFORM("All batches exhausted - current solution is the best!");
        return ompl::PlannerStatus::EXACT_SOLUTION;
      }
      else {
        /// No solution found
        OMPL_INFORM("All batches exhausted - no solution found!");
        return ompl::PlannerStatus::ABORT;
      }
    }

    if(ptc == true){
      OMPL_INFORM("Planner termination condition satisfied!");
      return ompl::PlannerStatus::TIMEOUT;
    }

    /// If it is the start of a new batch, get the next Batch
    /// And reset current inflation factor
    if(mIsInitSearchBatch) {
      mCurrentAlpha = 0.0;
      
      std::function<bool(Vertex)> pruneFunction = 
        std::bind(&BatchingPOMP::vertexPruneFunction, &BatchingPOMP,std::placeholders::_1);
      mBatchingPtr->nextBatch(pruneFunction, *mVertexNN)
    }

    /// Now do the search
    std::map<Vertex,Vertex> startPreds;
    std::map<Vertex,double> startDist;
    std::map<Vertex,double> startFValue;
    std::map<Vertex,boost::default_color_type> colorMap;

    try
    {
      // TODO - Make this more efficient!
      mNumSearches++;
      WeightMap wm;
      std::unique_ptr<boost::astar_heuristic<Graph, double>> heuristicFn;

      if(mIsInitSearchBatch) {
        /// Log of probability is cost - no heuristic
        wm = LogProbMap(*this);
        heuristicFn.reset(new zero_heuristic<Graph,double>());
      }
      else {
        /// Expected cost is cost
        wm = ExpWeightMap(*this);
        heuristicFn.reset(new zero_heuristic<Graph,double>(*this));
      }

      mIsInitSearchBatch = false; // Either way, make it false

      /// TODO : Find a less hacky way to do this
      if(mBatchingType == "single") {
        boost::astar_search(
          g,
          mStartVertex,
          *heuristicFn,
          throw_visitor(mGoalVertex),
          boost::make_assoc_property_map(startPreds),
          boost::make_assoc_property_map(startFValue),
          boost::make_assoc_property_map(startDist),
          wm,
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
          wm,
          boost::get(boost::vertex_index, g),
          boost::make_assoc_property_map(colorMap),
          std::less<double>(),
          boost::closed_plus<double>(std::numeric_limits<double>::max()),
          std::numeric_limits<double>::max(),
          double()
        );
      }
    }
    catch (const throw_visitor_exception & ex)
    {
    }


    if(startDist[mGoalVertex] == std::numeric_limits<double>::max()) {
      /// Did not find a path
      OMPL_INFORM("No further feasible paths to find in roadmap!");
      if(mBestPathCost < std::numeric_limits<double>::max()) {
        return ompl::PlannerStatus::EXACT_SOLUTION;
      }
      else {
        return ompl::PlannerStatus::ABORT;
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

    std::reverse(ePath.begin(), ePath.end());

    /// If path is free, set current cost to it, increment alpha
    bool pathBlocked{checkAndUpdatePathBlocked(ePath)};

    if(!pathBlocked) {
      mIsPathfound = true;
      currSolnCost = getPathDistance(ePath);

      if(mCurrentAlpha >= 1.0) {
        mIsInitSearchBatch = true;
      }
      else {
        mCurrentAlpha = std::min(mCurrentAlpha+mIncrement, 1.0);
      }
    }

    /// Path has been checked either way so update
    updateAffectedEdgeWeights();
  }


  /// Has found an improved feasible solution -> now report
  double improvementFactor{mBestPathCost/currSolnCost};
  mBestPathCost = currSolnCost;
  mCurrBestPath = ePath;

  ompl::geometric::PathGeometric* path = new ompl::geometric::PathGeometric(si_);
  VertexIndexMap vertex_id = get(boost::vertex_index,g);

  path->append(g[mStartVertex].v_state->state);

  for(auto e : mCurrBestPath)
  {
    path->append(g[target(e,g)].v_state->state);
  }

  pdef_->addSolutionPath(ompl::base::PathPtr(path));

  /// Prune vertices using current solution cost if sufficiently improved
  if(improvementFactor > BatchingPOMP::PRUNETHRESHOLD) {
    std::function<bool(Vertex)> admissibleFunction = 
        std::bind(&BatchingPOMP::isVertexAdmissible, &BatchingPOMP,std::placeholders::_1);
    mBatchingPtr->pruneVertices(admissibleFunction);
  }

  return ompl::base::PlannerStatus::APPROXIMATE_SOLUTION;

}


} // namespace batching_pomp