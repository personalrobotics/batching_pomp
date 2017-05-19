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

#ifndef BATCHING_POMP_BATCHING_MANAGER_HPP_
#define BATCHING_POMP_BATCHING_MANAGER_HPP_

#include <functional>
#include <ompl/base/StateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>
#include <boost/graph/adjacency_list.hpp>
#include "batching_pomp/util/RoadmapFromFile.hpp"

namespace batching_pomp {
namespace batching {

/// Abstract class that represents the batching strategy used for the
/// planning algorithm. 
template<class Graph, class VStateMap, class StateCon, class EDistance=StateCon>
class BatchingManager
{

typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::vertex_iterator VertexIter;
typedef typename GraphTypes::vertex_descriptor Vertex;

public:

  BatchingManager(const ompl::base::StateSpacePtr _space,
                  VStateMap _stateMap,
                  std::string _roadmapFileName,
                  Graph& _currentRoadmap)
  : mCurrentRoadmap{_currentRoadmap}
  , mNumBatches{0u}
  , mExhausted{false}
  , mCurrRadius{0.0}
  {
    auto file_roadmap_ptr = std::make_shared<
      batching_pomp::util::RoadmapFromFile<Graph,VStateMap,StateCon,EDistance>>
      (_space,_roadmapFileName);
    file_roadmap_ptr->generate(mFullRoadmap,_stateMap);
    mNumVertices = num_vertices(mFullRoadmap);
  }

  virtual ~BatchingManager() = default;

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
  
  void pruneVertices(std::function<bool(Vertex)>& _pruneFunction)
  {
    /// TODO: Check if this preserves iterator stability
    VertexIter next,vi,vi_end;
    boost::tie(vi,vi_end) = vertices(mCurrentRoadmap);

    for(next=vi; vi!=vi_end; vi=next)
    {
      if(_pruneFunction(*vi)) {
        remove_vertex(*vi,mCurrentRoadmap);
      }
    }
  }

  /// Nearest neighbour member may be nullptr (for single batch case)
  virtual void nextBatch(std::function<bool(Vertex)>& _pruneFunction,
                         ompl::NearestNeighbors<Vertex>& _vertexNN) = 0;


protected:

  Graph mFullRoadmap;
  Graph & mCurrentRoadmap;

  unsigned int mNumBatches;
  unsigned int mNumVertices;
  bool mExhausted;
  double mCurrRadius;

};

} // namespace batching
} // namespace batching_pomp

# endif //BATCHING_POMP_BATCHING_MANAGER_HPP_