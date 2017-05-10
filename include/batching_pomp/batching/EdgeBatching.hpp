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
               std::function<double(unsigned int)>& _initRadiusFn,
               double _radiusInflFactor)
  : BatchingManager(_space,_stateMap,_roadmapFileName)
  , mCurrRadius{_initRadiusFn(mNumVertices)}
  , mRadiusInflFactor{_radiusInflFactor}
  {
  }


private:

  double mRadiusInflFactor;
  double mCurrRadius;



};
} //namespace batching
} //namespace batching_pomp

#endif //BATCHING_POMP_EDGE_BATCHING_HPP_