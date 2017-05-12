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

template<class Graph, class VStateMap, class StateCon>
class EdgeBatching : public virtual BatchingManager<Graph, VStateMap, StateCon>
{

typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::vertex_iterator VertexIter;
typedef typename GraphTypes::vertex_descriptor Vertex;

public:

  EdgeBatching(const ompl::base::StateSpacePtr _space,
               VStateMap _stateMap,
               std::string _roadmapFileName,
               Graph& _currentRoadmap,
               double _radiusInflFactor,
               std::function<double(unsigned int)> _initRadiusFn,
               double _maxRadius
               )
  : BatchingManager(_space,_stateMap,_roadmapFileName,_currentRoadmap)
  , mRadiusInflFactor{_radiusInflFactor}
  , mInitRadius{_initRadiusFn(mNumVertices)}
  , mCurrRadius{mInitRadius}
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

  //////////////////////////////////////////////////
  /// Overriden methods
  void nextBatch(const std::function<bool(Vertex)>& _isAdmissible,
                 ompl::NearestNeighbors<Vertex>* _vertexNN) override
  {

    if(mExhausted){
      OMPL_INFORM("Batching exhausted! No updates with nextBatch!");
      return;
    }

    if(!_vertexNN){
      throw std::runtime_error("Nearest neighbour structure pointer is empty!")
    }

    OMPL_INFORM("New Edge Batch called!");
    ++mNumBatches;

    if(mNumBatches == 1u)
    {
      VertexIter vi, vi_end;
      std::vector<Vertex> vertex_vector(mNumVertices);
      size_t idx = 0;

      /// TODO : This assumes ALL vertices will be admissible (vertex_vector full)
      for(boost::tie(vi,vi_end)=vertices(mFullRoadmap); vi!=vi_end; ++vi)
      {
        if(_isAdmissible(mFullRoadmap[*mCurrVertex])) {
          Vertex newVertex = boost::add_vertex(mCurrentRoadmap);
          mCurrentRoadmap[newVertex].v_state = mFullRoadmap[*vi].v_state;
          vertex_vector[idx++] = newVertex;
        }
      }
      if(idx != mNumVertices){
        throw std::runtime_error("Not all vertices were admissible for edge batching!");
      }

      _vertexNN->add(vertex_vector);  
    }
    else {
      mCurrRadius = mCurrRadius * mRadiusInflFactor;
    }

    if(mCurrRadius > mMaxRadius)
    {
      mCurrRadius = mMaxRadius;
      mExhausted = true;
    }

  }

 
private:

  double mRadiusInflFactor;
  double mInitRadius;
  double mCurrRadius;
  double mMaxRadius;


};
} //namespace batching
} //namespace batching_pomp

#endif //BATCHING_POMP_EDGE_BATCHING_HPP_