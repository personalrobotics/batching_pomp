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

#ifndef BATCHING_POMP_BATCHINGPOMP_HPP_
#define BATCHING_POMP_BATCHINGPOMP_HPP_

#include <exception>
#include <memory>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

#include <ompl/base/Planner.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>

#include "batching.hpp"
#include "cspacebelief.hpp"
#include "util.hpp"


namespace batching_pomp {

/// Composite struct to associate the state space with each state
struct StateCon
{
  const ompl::base::StateSpacePtr space;
  ompl::base::State * state;
  StateCon(const ompl::base::StateSpacePtr _space)
  : space{_space}
  , state{space->allocState()} {}
  ~StateCon() {space->freeState(this->state); }
};

typedef std::shared_ptr<StateCon> StateConPtr;


class BatchingPOMP : public ompl::base::Planner
{

public:

  static const int FREE{1};
  static const int BLOCKED{-1};
  static const int UNKNOWN{0};
  static constexpr double DEFAULTPRUNETHRESHOLD = 1.1; //Threshold of old-cost/new-cost beyond which pruning should be done

  struct VProps
  {
    StateConPtr v_state;
  };

  struct EProps
  {
    double distance;
    double probFree;
    int blockedStatus;
  };

  /// Graph definitions
  typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, VProps, EProps> Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::vertex_iterator VertexIter;
  typedef boost::graph_traits<Graph>::edge_descriptor Edge;
  typedef boost::graph_traits<Graph>::edge_iterator EdgeIter;
  typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;

  /// Boost property maps needed by planner
  typedef boost::property_map<Graph, StateConPtr VProps::*>::type VPStateMap;
  typedef boost::property_map<Graph, double EProps::*>::type EPDistanceMap;
  typedef boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;
  //typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

  /// Public variables for the planner
  const ompl::base::StateSpacePtr mSpace;

  std::shared_ptr< batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > mBatchingPtr;

  std::shared_ptr< cspacebelief::Model<cspacebelief::BeliefPoint> > mBeliefModel;

  Graph g;

  /// OMPL methods
  BatchingPOMP(const ompl::base::SpaceInformationPtr & si);

  /// Constructor for non-default stuff
  BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
               std::shared_ptr< batching::BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > _batchingPtr,
               std::shared_ptr< cspacebelief::Model<cspacebelief::BeliefPoint> > _beliefModel,
               std::unique_ptr< util::Selector<Graph> > _selector,
               double _increment,
               double _startGoalRadius,
               double _pruneThreshold,
               const std::string& _roadmapFileName);

  ~BatchingPOMP(void);

  /// Setters and Getters
  double getCurrentAlpha() const;
  double getCurrentBestCost() const;
  Vertex getStartVertex() const;
  Vertex getGoalVertex() const;
  bool isInitSearchBatch() const;
  double getIncrement() const;
  void setIncrement(double _decrement);
  double getStartGoalRadius() const;
  void setStartGoalRadius(double _startGoalRadius);
  double getPruneThreshold() const;
  void setPruneThreshold(double _pruneThreshold);
  std::string getGraphType() const;
  void setGraphType(const std::string& _graphType);
  std::string getBatchingType() const;
  void setBatchingType(const std::string& _selectorType);
  std::string getSelectorType() const;
  void setSelectorType(const std::string& _selectorType);
  std::string getRoadmapFileName() const;
  void setRoadmapFileName(const std::string& _roadmapFileName);

  /// Evaluation methods
  inline unsigned int getNumEdgeChecks(){ return mNumEdgeChecks;}
  inline unsigned int getNumCollChecks(){ return mNumCollChecks;}
  inline unsigned int getNumSearches(){ return mNumSearches;}
  inline double getLookupTime(){ return mLookupTime;}

  /// OMPL required methods
  void setProblemDefinition(const ompl::base::ProblemDefinitionPtr & pdef);
  ompl::base::PlannerStatus solve(const ompl::base::PlannerTerminationCondition & ptc);

  /// Public helper methods
  double vertexDistFun(const Vertex& u, const Vertex& v) const;
  double computeAndSetEdgeFreeProbability(const Edge& e);
  bool checkAndSetEdgeBlocked(const Edge& e);

private:

  /// Planning helpers
  std::unique_ptr<util::Selector<Graph>> mSelector;
  std::unique_ptr<ompl::NearestNeighborsGNAT<Vertex>> mVertexNN;
  std::function<double(unsigned int)> mRadiusFun;
  util::BisectPerm mBisectPermObj;
  Vertex mStartVertex;
  Vertex mGoalVertex;
  std::vector<Edge> mEdgesToUpdate;
  std::vector<Edge> mCurrBestPath;

  /// Planner parameters
  double mCurrentAlpha;
  double mIncrement;
  double mStartGoalRadius;
  double mPruneThreshold;
  bool mIsInitSearchBatch;
  bool mIsPathFound;
  double mCheckRadius;
  double mBestPathCost;
  std::string mGraphType; // Optional - to be used for non-CPP level calls
  std::string mBatchingType; // Optional - to be used for non-CPP level calls
  std::string mSelectorType; // Optional - to be used for non-CPP level calls
  std::string mRoadmapName;

  /// For planner evaluation
  unsigned int mNumEdgeChecks;
  unsigned int mNumCollChecks;
  unsigned int mNumSearches;
  double mLookupTime;

  /// Private helper methods
  double haltonRadiusFun(unsigned int n) const;
  double rggRadiusFun(unsigned int n) const;
  void updateAffectedEdgeWeights();
  void addAffectedEdges(const Edge& e);
  bool checkAndUpdatePathBlocked(const std::vector<Edge>& _ePath);
  bool isVertexInadmissible(const Vertex& v) const ;
  bool vertexPruneFunction(const Vertex& v) const;
  double getPathDistance(const std::vector<Edge>& _ePath) const;
};

} // namespace batching_pomp
#endif //BATCHING_POMP_BATCHINGPOMP_HPP_