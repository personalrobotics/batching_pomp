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
#include "batching_pomp/include/BatchingPOMP.hpp"

namespace batching_pomp{

using batching_pomp::cspacebelief::BeliefPoint;

/// Default constructor for CPP-level usage
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
                           std::unique_ptr<BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>> _batchingPtr,
                           std::unique_ptr< Model<BeliefPoint> _beliefModel,
                           std::unique_ptr<Selector<Graph>> _selector,
                           double _searchInflFactor,
                           double _decrement,
                           double _startGoalRadius,
                           const std::string& _graphType,
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
, mGraphType{_graphType}
, mBatchingType{""}
, mRoadmapName{_roadmapFileName}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mLookupTime{0.0}
{
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
  mBeliefModel = std::make_unique<batching_pomp::cspacebelief::KNNModel>(15,2.0,0.5,0.25,bpDistFun);




}


BatchingPOMP::~BatchingPOMP()
{
}

////////////////////////////////////////////////////////////////////
/// Setters and Getters
double getSearchInflationFactor() const
{
  return mSearchInflationFactor;
}

void setSearchInflationFactor(double _searchInflFactor)
{
  mSearchInflationFactor = _searchInflFactor;
}

double getDecrement() const
{
  return mDecrement;
}

void setDecrement(double _decrement)
{
  mDecrement = _decrement;
}

double getStartGoalRadius() const
{
  return mStartGoalRadius;
}

void setStartGoalRadius(double _startGoalRadius)
{
  mStartGoalRadius = _startGoalRadius;
}

double getCheckRadius() const
{
  return mCheckRadius;
}

void setCheckRadius(double _checkRadius)
{
  mCheckRadius = _checkRadius;
}

std::string getGraphType() const
{
  return mGraphType;
}

void setGraphType(const std::string& _graphType)
{
  mGraphType = _graphType;
}

std::string getBatchingType() const
{
  return mBatchingType;
}

void setBatchingType(const std::string& _batchingType)
{
  mBatchingType = _batchingType;
}

std::string getRoadmapFileName() const
{
  return mRoadmapName;
}

void setRoadmapFileName(const std::string& _roadmapFileName)
{
  mRoadmapName = _roadmapFileName;
}







} // namespace batching_pomp