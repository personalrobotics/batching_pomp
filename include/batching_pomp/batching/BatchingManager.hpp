/***********************************************************************
Copyright (c) 2017, Shushman Choudhury
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

#ifndef BATCHING_POMP_BATCHING_MANAGER_HPP_
#define BATCHING_POMP_BATCHING_MANAGER_HPP_

#include <functional>
#include <ompl/base/StateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include <boost/graph/adjacency_list.hpp>
#include "batching_pomp/util/RoadmapFromFile.hpp"

namespace batching_pomp {
namespace batching {


//! Abstract class that represents the batching strategy used for the planning algorithm

/// The batching manager maintains responsibility
/// for the complete roadmap and periodically adds batches to the 
/// roadmap of the planner. This decoupling allows the complete roadmap
/// to be loaded separately before the planning query arrives.
/// @tparam Graph The type of boost graph used for the roadmaps
/// @tparam VStateMap The type of boost property map for vertex states
/// @tparam StateCon The wrapper type for an ompl state
/// @tparam EDistance The type of property map for edge lengths
template<class Graph, class VStateMap, class StateCon, class EDistance>
class BatchingManager
{

typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::vertex_iterator VertexIter;
typedef typename GraphTypes::vertex_descriptor Vertex;

public:

  /// \param[in] _space The state space of the planner
  /// \param[in] _stateMap The property map for states of vertices
  /// \param[in] _roadmapFileName The full path to roadmap .graphml file
  /// \param[in] _fullRoadmap The complete roadmap loaded by the batching manager
  /// \param[in] _currentRoadmap The roadmap that the planner is currently searching
  BatchingManager(const ompl::base::StateSpacePtr _space,
                  VStateMap _stateMap,
                  std::string _roadmapFileName,
                  Graph& _fullRoadmap,
                  Graph& _currentRoadmap)
  : mFullRoadmap(_fullRoadmap)
  , mCurrentRoadmap(_currentRoadmap)
  , mNumBatches{0u}
  , mExhausted{false}
  , mCurrRadius{0.0}
  {
    auto file_roadmap_ptr = std::make_shared<
      batching_pomp::util::RoadmapFromFile<Graph,VStateMap,StateCon,EDistance>>
      (_space,_roadmapFileName);
    file_roadmap_ptr->generateVertices(mFullRoadmap,_stateMap);
    mNumVertices = num_vertices(mFullRoadmap);
  }

  /// Corresponding constructor for single batching strategy
  /// \param[in] _space The state space of the planner
  /// \param[in] _stateMap The property map for states of vertices
  /// \param[in] _distanceMap The property map for the distance of edges
  /// \param[in] _roadmapFileName The full path to roadmap .graphml file
  /// \param[in] _fullRoadmap The complete roadmap loaded by the batching manager
  /// \param[in] _currentRoadmap The roadmap that the planner is currently searching
  BatchingManager(const ompl::base::StateSpacePtr _space,
                  VStateMap _stateMap,
                  EDistance _distanceMap,
                  std::string _roadmapFileName,
                  Graph& _fullRoadmap,
                  Graph& _currentRoadmap)
  : mFullRoadmap(_fullRoadmap)
  , mCurrentRoadmap(_currentRoadmap)
  , mNumBatches{0u}
  , mExhausted{false}
  , mCurrRadius{0.0}
  {
    auto file_roadmap_ptr = std::make_shared<
      batching_pomp::util::RoadmapFromFile<Graph,VStateMap,StateCon,EDistance>>
      (_space,_roadmapFileName);
    file_roadmap_ptr->generateVertices(mFullRoadmap,_stateMap);
    file_roadmap_ptr->generateEdges(mFullRoadmap,_stateMap,_distanceMap);
    mNumVertices = num_vertices(mFullRoadmap);
  }

  virtual ~BatchingManager() = default;

  //////////////////////////////////////////////////
  // Getters
  unsigned int getNumBatches() const
  {
    return mNumBatches;
  }

  unsigned int getNumVertices() const
  {
    return mNumVertices;
  }

  bool isExhausted() const
  {
    return mExhausted;
  }

  double getCurrentRadius() const
  {
    return mCurrRadius;
  }

  const ompl::base::State* getVertexState(const Vertex& v) const
  {
    return mFullRoadmap[v].v_state->state;
  }
  
  /// Implements the removal of all vertices, from the current roadmap,
  /// that satisfy a given pruning condition. Vertices are also
  /// removed from the nearest neighbour structure that maintains them.
  /// \param[in] _pruneFunction The function pointer that has the pruning condition
  /// \param[in] _vertexNN The nearest neighbour structure for vertices
  void pruneVertices(const std::function<bool(const ompl::base::State*)>& _pruneFunction,
                     ompl::NearestNeighbors<Vertex>& _vertexNN)
  {
    VertexIter vi,vi_end;
    std::vector<Vertex> verticesToRemove;
    verticesToRemove.reserve(num_vertices(mCurrentRoadmap));
    size_t vRemoved{0};

    for(boost::tie(vi,vi_end) = vertices(mCurrentRoadmap); vi!=vi_end; ++vi)
    { 
      if(_pruneFunction(mCurrentRoadmap[*vi].v_state->state)) {
        verticesToRemove[vRemoved++] = *vi;
      }
    }

    // Now remove vertices lined up
    for(size_t i=0; i<vRemoved; i++)
    {
      _vertexNN.remove(verticesToRemove[i]);
      clear_vertex(verticesToRemove[i],mCurrentRoadmap);
    }
  }

  /// Make any updates to batching manager with newest solution cost
  /// \param[in] _newSolnCost The current best solution cost with which
  ///                         to update the batching manager parameters
  virtual void updateWithNewSolutionCost(double _newSolnCost) = 0;

  /// Generate the next batch and update the current roadmap with it
  /// \param[in] _pruneFunction The rejection sampling checker for new samples
  /// \param[in] _vertexNN The nearest neighbour manager for roadmap vertices
  ///                       that is updated with the latest batch of samples
  virtual void nextBatch(const std::function<bool(const ompl::base::State*)>& _pruneFunction,
                         ompl::NearestNeighbors<Vertex>& _vertexNN) = 0;


protected:

  /// The entire roadmap that will be used for planning
  Graph & mFullRoadmap;

  /// The roadmap currently being searched by the planner
  Graph & mCurrentRoadmap;

  /// The number of batches added till now
  unsigned int mNumBatches;

  /// The number of vertices in the entire roadmap
  unsigned int mNumVertices;

  /// The flag that maintains whether the complete roadmap has been fully searched or not.
  bool mExhausted;

  /// The radius of connectivity for the current batch
  double mCurrRadius;

};

} // namespace batching
} // namespace batching_pomp

# endif //BATCHING_POMP_BATCHING_MANAGER_HPP_