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
#include <cmath>
#include "batching_pomp/include/BatchingPOMP.hpp"

namespace batching_pomp{

using batching_pomp::cspacebelief::BeliefPoint;

/// Default constructor for CPP-level usage
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
                           std::unique_ptr< BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > _batchingPtr,
                           std::unique_ptr< Model<BeliefPoint> > _beliefModel,
                           std::unique_ptr< Selector<Graph> > _selector,
                           double _searchInflFactor,
                           double _decrement,
                           double _startGoalRadius,
                           const std::string& _roadmapFileName)
: ompl::base::Planner(si,"BatchingPOMP")
, mSpace(si->getStateSpace()),
, mBatchingPtr{std::move(_batchingPtr)}
, mBeliefModel{std::move(_beliefModel)}
, mSelector{std::move(_selector)}
, mSearchInflationFactor{_searchInflFactor}
, mDecrement{_decrement}
, mStartGoalRadius{_startGoalRadius}
, mIsInitSearchBatch{false}
, mIsPathfound{false}
, mCheckRadius{0.5*space->getLongestValidSegmentLength()}
, mBestCost{std::numeric_limits<double>::max()}
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
, mSpace(si->getStateSpace()),
, mSearchInflationFactor{10.0}
, mDecrement{0.0}
, mStartGoalRadius{0.0}
, mIsInitSearchBatch{false}
, mIsPathfound{false}
, mCheckRadius{0.5*space->getLongestValidSegmentLength()}
, mBestCost{std::numeric_limits<double>::max()}
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
  Planner::declareParam<double>("inflation", this, &BatchingPOMP::setSearchInflationFactor, &BatchingPOMP::getSearchInflationFactor);
  Planner::declareParam<double>("decrement", this, &BatchingPOMP::setDecrement, &BatchingPOMP::getDecrement);
  Planner::declareParam<double>("start_goal_radius", this, &BatchingPOMP::setStartGoalRadius, &BatchingPOMP::getStartGoalRadius);
  Planner::declareParam<std::string>("graph_type", this, &BatchingPOMP::setSearchInflationFactor, &BatchingPOMP::getSearchInflationFactor);
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
double BatchingPOMP::getSearchInflationFactor() const
{
  return mSearchInflationFactor;
}

void BatchingPOMP::setSearchInflationFactor(double _searchInflFactor)
{
  mSearchInflationFactor = _searchInflFactor;
}

double BatchingPOMP::getDecrement() const
{
  return mDecrement;
}

void BatchingPOMP::setDecrement(double _decrement)
{
  mDecrement = _decrement;
}

double BatchingPOMP::getStartGoalRadius() const
{
  return mStartGoalRadius;
}

void BatchingPOMP::setStartGoalRadius(double _startGoalRadius)
{
  mStartGoalRadius = _startGoalRadius;
}

double BatchingPOMP::getCheckRadius() const
{
  return mCheckRadius;
}

void BatchingPOMP::setCheckRadius(double _checkRadius)
{
  mCheckRadius = _checkRadius;
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


/// Private helper methods
double BatchingPOMP::vertexDistFun(const Vertex& u, const Vertex& v) const
{
  return mSpace->distance(g[u].v_state->state, g[v].v_state->state);
}

double BatchingPOMP::haltonRadiusFun(unsigned int n) const
{
  auto dimDbl = static_cast<double>(mSpace->getDimension());
  // TODO : Add the formula here
}

double BatchingPOMP::rggRadiusFun(unsigned int n) const
{
  auto dimDbl = static_cast<double>(mSpace->getDimension());
  // TODO : Add the formula here
}

//TODO : USE STATECONPTR HERE!
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
    ompl::base::State* tempState{mSpace->allocState()};
    mSpace->interpolate(startState, endState,
      1.0*(1+i)/(nStates+1), tempState);

    BeliefPoint query(tempState, -1.0);
    double collProb = mBeliefModel->estimate(query);
    result *= (1.0 - coll_prob);
  }

  g[e].probFree = result;
  return result;
}


bool checkAndSetEdgeBlocked(const Edge& e)
{
  /// March along edge states with highest resolution
  mNumEdgeChecks++;

  auto validityChecker = si_->getStateValidityChecker();
  
  auto startState = g[source(e,g)].v_state->state;
  auto endState = g[target(e,g)].v_state->state;
  unsigned int nStates{std::floor(g[e].distance / (2.0*mCheckRadius))};

  const std::vector< std::pair<int,int> > & order = bisectPermObj_.get(nStates);

  for(unsigned int i = 0; i < nStates; i++)
  {
    ompl::base::State* tempState{mSpace->allocState()};
    mSpace->interpolate(startState, endState,
      1.0*(1+order[ui].first)/(nStates+1), tempState);

    mNumCollChecks++;
    /// Check and add to belief model
    if(validityChecker->isValid(tempState) == false) {
      BeliefPoint toAdd(tempState,1.0);
      mBeliefModel->addPoint(toAdd);
      g[e].blockedStatus = BatchingPOMP::BLOCKED;
      return true;
    }
    else {
      BeliefPoint toAdd(tempState,0.0);
      mBeliefModel->addPoint(toAdd);
    }
  }

  g[e].blockedStatus = BatchingPOMP::FREE;
  return false;
}






} // namespace batching_pomp