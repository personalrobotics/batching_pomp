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


#include "batching.hpp"
#include "cspacebelief.hpp"
#include "util.hpp"


namespace batching_pomp {

/// Composite struct to associate the state space with each state
struct StateCon
{
  const ompl::base::StateSpacePtr space;
  ompl::base::State * state;
  StateCon(ompl::base::StateSpacePtr space)
  : space(space)
  , state(space->allocState()) {}
  ~StateCon() {space->freeState(this->state); }
};


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
    if(u == mVThrow)
      throw throw_visitor_exception();
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
  ompl::NearestNeighbors<BatchingPOMP::Vertex>& mVertexNN;
  double mCurrRadius;
  std::function<double(const BatchingPOMP::Vertex&, const BatchingPOMP::Vertex&)> mVertexDistFun;

  neighbours_visitor(ompl::NearestNeighbors<BatchingPOMP::Vertex>& _vertexNN, 
                     double _currRadius)
  : mVertexNN{_vertexNN}
  , mCurrRadius{_currRadius} 
  , mVertexDistFun{mVertexNN->getDistanceFunction()} {}
  inline void initialize_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void discover_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g) {}
  inline void examine_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph& g)
  {
    std::vector<BatchingPOMP::Vertex> vertexNbrs;
    mVertexNN.nearestR(u,mCurrRadius,vertexNbrs);

    // Now iterate through neighbors and check if edge exists between
    // u and neighbour. If it does not exist, create it AND set its
    // distance (one-time computation) based on provided distance function
    for(Vertex nbr : vertexNbrs){
      if(!edge(u,nbr,g).second){
        std::pair<Edge,bool> new_edge = add_edge(u,nbr,g);
        g[new_edge.first].distance = mVertexDistFun(source(new_edge,g), target(new_edge,g));
      }
    }
  }
  inline void examine_edge(BatchingPOMP::Edge e, const BatchingPOMP::Graph& g) {}
  inline void edge_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void edge_not_relaxed(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void black_target(BatchingPOMP::Edge e, const BatchingPOMP::Graph & g) {}
  inline void finish_vertex(BatchingPOMP::Vertex u, const BatchingPOMP::Graph & g) {}
};



} // namespace batching_pomp
#endif //BATCHING_POMP_BATCHINGPOMP_HPP_