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
#include "GraphTypes.hpp"

namespace batching_pomp {


/// The OMPL Planner class that implements the algorithm
class BatchingPOMP : public ompl::base::Planner
{

public:

  /// The pointer to the OMPL state space 
  const ompl::base::StateSpacePtr mSpace;

  /// The pointer to the batching manager instance to be used by the planner
  std::shared_ptr< batching::BatchingManager > mBatchingPtr;

  /// The pointer to the C-space belief model instance to be used by the planner
  std::shared_ptr< cspacebelief::Model<cspacebelief::BeliefPoint> > mBeliefModel;

  std::unique_ptr<ompl::NearestNeighborsGNAT<Vertex>> mVertexNN;

  /// The roadmap that will be continuously searched and updated.
  Graph* g;

  /// The large, dense roadmap that remains unchanged after being loaded once
  Graph* full_g;

  BatchingPOMP(const ompl::base::SpaceInformationPtr & si);

  /// \param[in] si The OMPL space information manager
  /// \param[in] _batchingPtr The pointer to the constructed batching manager that
  ///            the planner should use
  /// \param[in] _beliefModel The pointer to the constructed C-space belief model that
  ///            the planner should use
  /// \param[in] _selector The pointer to the constructed edge selector that 
  ///            the planner should use.
  /// \param[in] _roadmapFileName The path to the .graphml file that encodes the roadmap vertices
  /// \param[in] _startGoalRadius (Only for Single Batching) The radius to connect start
  ///            and goal vertices to the roadmap
  /// \param[in] _increment The increment in the value of alpha after each POMP search
  /// \param[in] _pruneThreshold The fractional change in cost above which pruning of samples should be done
  BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
               std::shared_ptr< batching::BatchingManager >& _batchingPtr,
               std::shared_ptr< cspacebelief::Model<cspacebelief::BeliefPoint> >& _beliefModel,
               std::shared_ptr< util::Selector<Graph>>& _selector,
               const std::string& _roadmapFileName,
               double _startGoalRadius,
               double _increment = 0.2,
               double _pruneThreshold = 0.05);

  ~BatchingPOMP(void);

  // Setters and Getters
  /// Current value of alpha being used by POMP
  double getCurrentAlpha() const;
  /// Cost of current best solution to goal
  double getCurrentBestCost() const;
  /// The ID of the start vertex
  Vertex getStartVertex() const;
  /// The ID of the goal vertex
  Vertex getGoalVertex() const;
  /// Whether the current POMP search is the first of that batch (alpha = 0)
  bool isInitSearchBatch() const;
  /// The value with which alpha is incremented
  double getIncrement() const;
  void setIncrement(double _decrement);
  double getStartGoalRadius() const;
  void setStartGoalRadius(double _startGoalRadius);
  double getPruneThreshold() const;
  void setPruneThreshold(double _pruneThreshold);
  /// The kind of sequence used to generate samples ("halton" or "rgg")
  std::string getGraphType() const;
  void setGraphType(const std::string& _graphType);
  std::string getBatchingType() const;
  void setBatchingType(const std::string& _selectorType);
  std::string getSelectorType() const;
  void setSelectorType(const std::string& _selectorType);
  std::string getRoadmapFileName() const;
  void setRoadmapFileName(const std::string& _roadmapFileName);

  // Internal evaluation methods
  /// Number of edges evaluated thus far
  inline unsigned int getNumEdgeChecks(){ return mNumEdgeChecks;}
  /// Number of calls to collision checker made thus far
  inline unsigned int getNumCollChecks(){ return mNumCollChecks;}
  /// Number of roadmap searches done thus far
  inline unsigned int getNumSearches(){ return mNumSearches;}
  /// Number of model lookups made thus far
  inline unsigned int getNumLookups(){ return mNumLookups;}
  /// Total time spent doing model lookups
  inline double getLookupTime(){ return mLookupTime;}
  /// Total time spent doing searches
  inline double getSearchTime(){ return mSearchTime;}
  /// Total time spent doing collision checks
  inline double getCollCheckTime(){return mCollCheckTime;}

  // OMPL required methods
  void setProblemDefinition(const ompl::base::ProblemDefinitionPtr & pdef);
  ompl::base::PlannerStatus solve(const ompl::base::PlannerTerminationCondition & ptc);
  void setup();

  /// The distance function between roadmap vertices to be used
  /// by the nearest neighbour manager for vertices. Typically
  /// returns the distance between the underlying states of the space.
  /// \param[in] u,v The end-point vertices of the edge
  /// \return The distance between vertices
  double vertexDistFun(const Vertex& u, const Vertex& v) const;

  /// Given a new edge, initialize the embedded configurations
  /// along the edge using the resolution of the underlying space.
  /// This is a separate method so that it is called only when
  /// a new edge is created, i.e. just-in-time.
  /// \param[in] e The edge ID to initialize with configurations
  void initializeEdgePoints(const Edge& e);

  /// Compute the collision measure of an edge and assign it
  /// to the underlying property of the edge. Also return the value
  /// to the caller.
  /// \param[in] e The edge ID
  /// \return The computed collision measure of e as per the current model
  double computeAndSetEdgeCollisionMeasure(const Edge& e);

  /// Compute the collision measure of a potential edge
  /// WITHOUT states initialized underneath
  /// \param[in] e The edge ID
  /// \return The computed collision measure of e as per the current model
  double computeEdgeCollisionMeasureNoStates(const Vertex& u, const Vertex& v);

  /// Evaluate an edge to determine its collision status
  /// and assign it to the underlying property of the edge.
  /// \param[in] The edge ID to check for
  /// \return True or False depending on if the edge is in collision or not
  bool checkAndSetEdgeBlocked(const Edge& e);


  void assignBatchingRoadmap();

private:

  // Planning helpers
  std::shared_ptr< util::Selector<Graph>> mSelector;
  std::function<double(unsigned int)> mRadiusFun;
  util::BisectPerm mBisectPermObj;
  Vertex mStartVertex;
  Vertex mGoalVertex;
  std::vector<Edge> mCurrBestPath;
  std::set<Edge> mEdgesToUpdate;

  // Planner parameters
  double mCurrentAlpha;
  double mIncrement;
  double mStartGoalRadius;
  double mPruneThreshold;
  bool mIsInitSearchBatch;
  bool mIsPathFound;
  double mCheckRadius;
  double mBestPathCost;
  bool mUsingSelector;
  std::string mGraphType; // Optional - to be used for non-CPP level calls
  std::string mBatchingType; // Optional - to be used for non-CPP level calls
  std::string mSelectorType; // Optional - to be used for non-CPP level calls
  std::string mRoadmapName;

  // For planner evaluation
  unsigned int mNumEdgeChecks;
  unsigned int mNumCollChecks;
  unsigned int mNumSearches;
  unsigned int mNumLookups;
  double mLookupTime;
  double mCollCheckTime;
  double mSearchTime;


  // Private helper methods
  double haltonRadiusFun(unsigned int n) const;
  double rggRadiusFun(unsigned int n) const;
  bool checkAndUpdatePathBlocked(const std::vector<Edge>& _ePath);
  void addAffectedEdges(const Edge& e);
  void updateAffectedEdgeWeights();
  bool isVertexInadmissible(const ompl::base::State* vState) const ;
  bool vertexPruneFunction(const ompl::base::State* vState) const;
  double getPathDistance(const std::vector<Edge>& _ePath) const;
};

} // namespace batching_pomp
#endif //BATCHING_POMP_BATCHINGPOMP_HPP_