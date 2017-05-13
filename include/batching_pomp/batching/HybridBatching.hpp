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

#ifndef BATCHING_POMP_HYBRID_BATCHING_HPP_
#define BATCHING_POMP_HYBRID_BATCHING_HPP_

#include <exception>
#include <ompl/base/StateSpace.h>
#include <ompl/util/Console.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include "batching_pomp/batching/BatchingManager.hpp"

namespace batching_pomp {
namespace batching {

template<class Graph, class VStateMap, class StateCon>
class HybridBatching : public BatchingManager<Graph, VStateMap, StateCon> 
{

typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::vertex_iterator VertexIter;
typedef typename GraphTypes::vertex_descriptor Vertex;

using BatchingManager<Graph, VStateMap, StateCon>::mFullRoadmap;
using BatchingManager<Graph, VStateMap, StateCon>::mCurrentRoadmap;
using BatchingManager<Graph, VStateMap, StateCon>::mNumBatches;
using BatchingManager<Graph, VStateMap, StateCon>::mNumVertices;
using BatchingManager<Graph, VStateMap, StateCon>::mExhausted;

public:

  HybridBatching(const ompl::base::StateSpacePtr _space,
                 VStateMap _stateMap,
                 std::string _roadmapFileName,
                 Graph& _currentRoadmap,
                 unsigned int _initNumVertices,
                 double _vertInflFactor,
                 double _radiusInflFactor,
                 std::function<double(unsigned int)> _initRadiusFn,
                 double _maxRadius
                 )
  : BatchingManager<Graph, VStateMap, StateCon>(_space,_stateMap,_roadmapFileName,_currentRoadmap)
  , mNumVerticesAdded{0u}
  , mNextVertexTarget{_initNumVertices}
  , mVertInflFactor{_vertInflFactor}
  , mRadiusInflFactor{_radiusInflFactor}
  , mInitRadius{_initRadiusFn(_initNumVertices)}
  , mCurrRadius{mInitRadius}
  , mMaxRadius{_maxRadius}
  , mEdgeBatchingMode{false}
  , mRadiusFn{std::move(_initRadiusFn)}
  {
    boost::tie(mCurrVertex,mLastVertex) = vertices(mFullRoadmap);
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

  void setCurrentRadius(double _currentRadius)
  {
    mCurrRadius = _currentRadius;
  }

  double getCurrentRadius() const
  {
    return mCurrRadius;
  }

  void setMaxRadius(double _maxRadius)
  {
    mMaxRadius = _maxRadius;
  }

  double getMaxRadius() const
  {
    return mMaxRadius;
  }

  bool isInEdgeBatchingMode() const
  {
    return mEdgeBatchingMode;
  }

  //////////////////////////////////////////////////
  /// Overriden methods
  void nextBatch(const std::function<bool(Vertex)>& _pruneFunction,
                 ompl::NearestNeighbors<Vertex>* _vertexNN) override
  {

    if(mExhausted){
      OMPL_INFORM("Batching exhausted! No updates with nextBatch!");
      return;
    }

    if(!_vertexNN){
      throw std::runtime_error("Nearest neighbour structure pointer is empty!");
    }

    OMPL_INFORM("New Hybrid Batch called!");
    ++mNumBatches;

    if(!mEdgeBatchingMode)
    {
      std::vector<Vertex> vertex_vector;

      while(mNumVerticesAdded < mNextVertexTarget)
      {

        if(_pruneFunction(mFullRoadmap[*mCurrVertex])) {
          Vertex newVertex = boost::add_vertex(mCurrentRoadmap);
          mCurrentRoadmap[newVertex].v_state = mFullRoadmap[*mCurrVertex].v_state;
          vertex_vector.push_back(newVertex);
        }

        ++mCurrVertex;
        ++mNumVerticesAdded;

        /// If all samples added, next round will be edge batching mode
        if(mCurrVertex == mLastVertex) {
          mEdgeBatchingMode = true;
          break;
        }
      }

      if(vertex_vector.size()>0) {
        _vertexNN->add(vertex_vector);
      }
      
      mNextVertexTarget = static_cast<unsigned int>(mNextVertexTarget * mVertInflFactor);

      mCurrRadius = mRadiusFn(mNextVertexTarget);
    }
    else
    {
      mCurrRadius = mCurrRadius * mRadiusInflFactor;

      if(mCurrRadius > mMaxRadius)
      {
        mCurrRadius = mMaxRadius;
        mExhausted = true;
      }
    }

  }


private:

  unsigned int mNumVerticesAdded;
  unsigned int mNextVertexTarget;
  double mVertInflFactor;

  double mRadiusInflFactor;
  double mInitRadius;
  double mCurrRadius;
  double mMaxRadius;

  VertexIter mCurrVertex;
  VertexIter mLastVertex;

  bool mEdgeBatchingMode;
  std::function<double(unsigned int)> mRadiusFn;

};

} // namespace batching_pomp
} // namespace batching
#endif //BATCHING_POMP_HYBRID_BATCHING_HPP_