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

#include <ompl/base/Planner.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>

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
  : space(_space)
  , state(space->allocState()) {}
  ~StateCon() {space->freeState(this->state); }
};

typedef std::shared_ptr<StateCon> StateConPtr;


class BatchingPOMP : public ompl::base::Planner
{

  using batching_pomp::batching::BatchingManager;
  using batching_pomp::cspacebelief::BeliefPoint;
  using batching_pomp::cspacebelief::Model;
  using batching_pomp::cspacebelief::Selector;
  using batching_pomp::cspacebelief::BisectPerm;

  const static int FREE{1};
  const static int BLOCKED{-1};
  const static int UNKNOWN{0};
  const static double PRUNETHRESOLD{1.1}; //Threshold of old-cost/new-cost beyond which pruning should be done

public:

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
  typedef boost::property_map<Graph, boost::edge_weight_t>::type WeightMap;

  /// Public variables for the planner
  const ompl::base::StateSpacePtr mSpace;

  std::unique_ptr< BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > mBatchingPtr;

  std::unique_ptr< Model<BeliefPoint> > mBeliefModel;

  Graph g;

  /// OMPL methods
  BatchingPOMP(const ompl::base::SpaceInformationPtr & si);

  /// Constructor for non-default stuff
  BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
               std::unique_ptr< BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap> > _batchingPtr,
               std::unique_ptr< Model<BeliefPoint> _beliefModel,
               std::unique_ptr< Selector<Graph> > _selector,
               double _searchInflFactor,
               double _decrement,
               double _startGoalRadius,
               double _checkRadius,
               const std::string& _roadmapFileName);

  ~BatchingPOMP(void);

  /// Setters and Getters
  double getCurrentAlpha() const;
  double getCurrentBestCost() const;
  Vertex getStartVertex() const;
  Vertex getGoalVertex() const;
  bool isInitSearchBatch();
  void setSearchInflationFactor(double _searchInflFactor);
  double getIncrement() const;
  void setIncrement(double _decrement);
  double getStartGoalRadius() const;
  void setStartGoalRadius(double _startGoalRadius);
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
  double computeAndSetEdgeFreeProbability(const Edge& e);
  bool checkAndSetEdgeBlocked(const Edge& e);
  double vertexDistFun(const Vertex& u, const Vertex& v) const;

private:

  /// Planning helpers
  std::unique_ptr<Selector<Graph>> mSelector;
  std::unique_ptr<ompl::NearestNeighborsGNAT<BeliefPoint>> mVertexNN;
  std::function<double(unsigned int)> mRadiusFun;
  BisectPerm mBisectPermObj;
  Vertex mStartVertex;
  Vertex mGoalVertex;
  std::vector<Edge> mEdgesToUpdate;
  std::vector<Edge> mCurrBestPath;
  AStarVisitor mSearchVisitor;

  /// Planner parameters
  double mCurrentAlpha;
  double mIncrement;
  double mStartGoalRadius;
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
  bool isVertexAdmissible(const Vertex& v)
  bool vertexPruneFunction(const Vertex& v) const;
  double getPathDistance(const std::vector<Edge>& _ePath) const;
};




////////////////////////////////////////////////////////////////////
/// Helper classes

class throw_visitor_exception : public std::exception {};
class throw_visitor
{

public:
  BatchingPOMP::Vertex mVThrow;
  throw_visitor(BatchingPOMP::Vertex _vThrow): mVThrow{_vThrow} {}
  inline void initialize_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void discover_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void examine_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g)
  {
    if(u == mVThrow) {
      throw throw_visitor_exception();
    }
  }
  inline void examine_edge(BatchingPOMP::Edge e, const BatchingPOMP::Graph& g) {}
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
  : mPlanner{_planner}
  , mVertexNN{_vertexNN}
  , mCurrRadius{_currRadius}
  , mVThrow{_vThrow}
  , mVertexDistFun{mVertexNN.getDistanceFunction()} {}
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
    for(Vertex nbr : vertexNbrs){
      if(!edge(u,nbr,g).second){
        std::pair<Edge,bool> new_edge = add_edge(u,nbr,g);
        g[new_edge.first].distance = mVertexDistFun(source(new_edge,g), target(new_edge,g));
        g[new_edge.first].blockedStatus = BatchingPOMP::UNKNOWN;
        mPlanner.computeAndSetEdgeFreeProbability(new_edge.first);
      }
    }
  }
  inline void examine_edge(BatchingPOMP::Edge e, const BatchingPOMP::Graph& g) {}
  inline void edge_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void edge_not_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void black_target(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void finish_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph & g) {}
};

/// Weight map for expected edge cost
class ExpWeightMap
{
public:

  const BatchingPOMP& mPlanner;
  EdgeWeightMap(const BatchingPOMP& _planner)
  : mPlanner{_planner} {}
};

const double get(const EdgeWeightMap& _ewMap, const BatchingPOMP::Edge& e)
{
  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::BLOCKED) {
    return std::numeric_limits<double>::max();
  }

  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::FREE) {
    return _ewMap.mPlanner.g[e].distance;
  }

  double alpha{_ewMap.mPlanner->getCurrentAlpha()};
  double w_m{-std::log(_ewMap.mPlanner.g[e].probFree)};
  double w_l{dmap.pm_planner.g[e].distance};

  double exp_cost{alpha*w_l + (1-alpha)*w_m};

  return exp_cost;
}

/// Log probability weight map
class LogProbMap
{
public:

  const BatchingPOMP& mPlanner;
  LogProbMap(BatchingPOMP& _planner)
  : mPlanner{_planner} {}
};

const double get(const LogProbMap& _ewMap, const BatchingPOMP::Edge& e)
{
  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::BLOCKED) {
    return std::numeric_limits<double>::max();
  }

  if(_ewMap.mPlanner.g[e].blockedStatus == BatchingPOMP::FREE) {
    return 0.0;
  }

  return -std::log(g[e].probFree);
};

/// Euclidean distance heuristic for A-star search
template<class Graph, class CostType>
class exp_distance_heuristic : public boost::astar_heuristic<Graph, CostType>
{
public:
  exp_distance_heuristic(const BatchingPOMP& _planner)
  : mPlanner{_planner}

  CostType operator()(batching_pomp::BatchingPOMP::Vertex u)
  {
    return mPlanner->getCurrentAlpha() * mPlanner.vertexDistFun(u, mPlanner->getGoalVertex());
  }

private:
  const BatchingPOMP& mPlanner;
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


} // namespace batching_pomp
#endif //BATCHING_POMP_BATCHINGPOMP_HPP_