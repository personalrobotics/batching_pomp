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

#ifndef BATCHING_POMP_VERTEX_BATCHING_HPP_
#define BATCHING_POMP_VERTEX_BATCHING_HPP_

#include <exception>
#include <ompl/base/StateSpace.h>
#include <ompl/util/Console.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include "batching_pomp/batching/BatchingManager.hpp"

namespace batching_pomp {
namespace batching {

//! Derived class of BatchingManager that implements Vertex Batching.

/// Implements Vertex Batching, where batches of vertices are added
/// and the complete subgraph is searched over each time.
class VertexBatching : public BatchingManager 
{

public:

  /// \param[in] _initNumVertices The initial number of vertices to begin with
  /// \param[in] _vertInflFactor The factor by which to increase the number of
  ///                            vertices in each batch.
  VertexBatching(const ompl::base::StateSpacePtr _space,
                 std::string _roadmapFileName,
                 unsigned int _initNumVertices,
                 double _vertInflFactor
                 )
  : BatchingManager(_space,_roadmapFileName)
  , mNumVerticesAdded{0u}
  , mNextVertexTarget{_initNumVertices}
  , mVertInflFactor{_vertInflFactor}
  {
    BatchingManager::mBatchingType = "vertex";
    BatchingManager::mCurrRadius = std::numeric_limits<double>::max();
    boost::tie(mCurrVertex,mLastVertex) = vertices(BatchingManager::mFullRoadmap);
  }

  //////////////////////////////////////////////////
  /// Setters and Getters
  void setVertexInflationFactor(double _vertInflFactor)
  {
    mVertInflFactor = _vertInflFactor;
  }

  double getVertexInflationFactor() const
  {
    return mVertInflFactor;
  }

  //////////////////////////////////////////////////
  /// Overriden methods
  void updateWithNewSolutionCost(double _newSolnCost) override
  {
  }


  void nextBatch(const std::function<bool(const ompl::base::State*)>& _pruneFunction,
                 ompl::NearestNeighbors<Vertex>& _vertexNN) override
  {

    if(BatchingManager::mExhausted){
      OMPL_INFORM("Batching exhausted! No updates with nextBatch!");
      return;
    }

    OMPL_INFORM("New Vertex Batch called!");
    ++BatchingManager::mNumBatches;

    std::vector<Vertex> vertex_vector(mNextVertexTarget - mNumVerticesAdded);
    size_t idx{0};

    while(mNumVerticesAdded < mNextVertexTarget)
    {

      // Only add if best cost through vertex better than current solution
      if(!_pruneFunction(BatchingManager::mFullRoadmap[*mCurrVertex].v_state->state)) {
        Vertex newVertex{boost::add_vertex(BatchingManager::mCurrentRoadmap)};
        BatchingManager::mCurrentRoadmap[newVertex].v_state = BatchingManager::mFullRoadmap[*mCurrVertex].v_state;
        vertex_vector[idx++] = newVertex;
      }

      // Increment stuff
      ++mCurrVertex;
      ++mNumVerticesAdded;

      if(mCurrVertex == mLastVertex) {
        BatchingManager::mExhausted = true;
        break;
      }
    }

    // Update nearest neighbour structure with vertices
    // TODO : Is this the most efficient way?
    if(idx > 0u) {
      vertex_vector.resize(idx);
      _vertexNN.add(vertex_vector);
    }

    // Update size of next subgraph
    mNextVertexTarget = static_cast<unsigned int>(mNextVertexTarget * mVertInflFactor);
  }

private:

  unsigned int mNumVerticesAdded;
  unsigned int mNextVertexTarget;
  double mVertInflFactor;

  VertexIter mCurrVertex;
  VertexIter mLastVertex;

};

} //namespace batching
} //namespace batching_pomp

#endif //BATCHING_POMP_VERTEX_BATCHING_HPP_
