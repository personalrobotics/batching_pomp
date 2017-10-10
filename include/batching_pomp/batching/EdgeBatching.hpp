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

#ifndef BATCHING_POMP_EDGE_BATCHING_HPP_
#define BATCHING_POMP_EDGE_BATCHING_HPP_

#include <exception>
#include <ompl/base/StateSpace.h>
#include <ompl/util/Console.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include "batching_pomp/batching/BatchingManager.hpp"

namespace batching_pomp {
namespace batching {

//! Derived class of BatchingManager that implements Edge Batching.

/// Implements Edge Batching, where all vertices are considered
/// and batches of edges are added based on an increasing
/// radius of connectivity.
class EdgeBatching : public BatchingManager
{

public:

  /// \param[in] _radiusInflFactor The value by which to increase the radius for each new batch
  /// \param[in] _initRadiusFn The function which generates the initial radius based on
  ///                          the number of vertices
  /// \param[in] _maxRadius The maximum meaningful radius for edge batching.
  ///                        It is a function of the space bounds and the current solution quality.
  EdgeBatching(const ompl::base::StateSpacePtr _space,
               std::string _roadmapFileName,
               double _radiusInflFactor,
               std::function<double(unsigned int)> _initRadiusFn,
               double _maxRadius
               )
  : BatchingManager(_space,_roadmapFileName)
  , mRadiusInflFactor{_radiusInflFactor}
  , mInitRadius{_initRadiusFn(BatchingManager::mNumVertices)}
  , mMaxRadius{_maxRadius}
  {
  }

  //////////////////////////////////////////////////
  /// Setters and Getters
  void setRadiusInflationFactor(unsigned int _radiusInflFactor)
  {
    mRadiusInflFactor = _radiusInflFactor;
  }

  double getRadiusInflationFactor() const
  {
    return mRadiusInflFactor;
  }

  void setInitRadius(double _initRadius)
  {
    mInitRadius = _initRadius;
  }

  double getInitRadius() const
  {
    return mInitRadius;
  }

  void setMaxRadius(double _maxRadius)
  {
    mMaxRadius = _maxRadius;
  }

  double getMaxRadius() const
  {
    return mMaxRadius;
  }

  //////////////////////////////////////////////////
  /// Overriden methods
  void updateWithNewSolutionCost(double _newSolnCost) override
  {
    mMaxRadius = std::min(mMaxRadius,_newSolnCost);
  }

  void nextBatch(const std::function<bool(const ompl::base::State*)>& _pruneFunction,
                 ompl::NearestNeighbors<Vertex>& _vertexNN) override
  {

    if(BatchingManager::mExhausted){
      OMPL_INFORM("Batching exhausted! No updates with nextBatch!");
      return;
    }

    OMPL_INFORM("New Edge Batch called!");
    ++BatchingManager::mNumBatches;

    if(BatchingManager::mNumBatches == 1u)
    {
      VertexIter vi, vi_end;
      std::vector<Vertex> vertex_vector(BatchingManager::mNumVertices);
      size_t idx{0};

      for(boost::tie(vi,vi_end)=vertices(BatchingManager::mFullRoadmap); vi!=vi_end; ++vi)
      {
        if(!_pruneFunction(BatchingManager::mFullRoadmap[*vi].v_state->state)) {
          Vertex newVertex = boost::add_vertex(BatchingManager::mCurrentRoadmap);
          BatchingManager::mCurrentRoadmap[newVertex].v_state = 
              BatchingManager::mFullRoadmap[*vi].v_state;
          vertex_vector[idx++] = newVertex;
        }
      }

      // Truncate to actual number of samples to add
      vertex_vector.resize(idx);
      _vertexNN.add(vertex_vector);

      BatchingManager::mCurrRadius = mInitRadius;  
    }
    else {
      BatchingManager::mCurrRadius = BatchingManager::mCurrRadius * mRadiusInflFactor;
    }

    if(BatchingManager::mCurrRadius > mMaxRadius)
    {
      BatchingManager::mCurrRadius = mMaxRadius;
      BatchingManager::mExhausted = true;
    }
  }

 
private:

  double mRadiusInflFactor;
  double mInitRadius;
  double mMaxRadius;


};
} //namespace batching
} //namespace batching_pomp

#endif //BATCHING_POMP_EDGE_BATCHING_HPP_