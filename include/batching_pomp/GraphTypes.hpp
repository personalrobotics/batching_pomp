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

#ifndef BATCHING_POMP_GRAPHTYPES_HPP_
#define BATCHING_POMP_GRAPHTYPES_HPP_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/dynamic_property_map.hpp>

namespace batching_pomp {

/// Composite struct to associate the state space with each state
struct StateCon
{
  const ompl::base::StateSpacePtr space;
  ompl::base::State * state;
  StateCon(const ompl::base::StateSpacePtr _space)
  : space{_space}
  , state{space->allocState()} {}
  ~StateCon() {space->freeState(this->state); }
};
typedef std::shared_ptr<StateCon> StateConPtr;

////////////////////////////////////////////////////////////////
// Graph Types

/// The edge is known to be collision-free
const int FREE{1};

/// The edge is known to be in collision
const int BLOCKED{-1};

/// The collision status of the edge is unknown
const int UNKNOWN{0};

/// Properties associated with each roadmap vertex
struct VProps
{
  /// The underlying state of the vertex
  StateConPtr v_state;
};

/// Properties associated with each roadmap edge
struct EProps
{
  /// The length of the edge using the space distance metric
  double distance;

  /// The probability of collision of the edge
  double collMeasure;

  /// The collision status of the edge (free, blocked or unknown)
  int blockedStatus;

  /// Have the embedded points been initialized
  bool hasPoints;

  /// The set of states embedded in the edge
  std::vector< StateConPtr > edgeStates;
};

// Graph definitions
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, VProps, EProps> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertex_iterator VertexIter;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIter;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;

// Boost property maps needed by planner
typedef boost::property_map<Graph, StateConPtr VProps::*>::type VPStateMap;
typedef boost::property_map<Graph, double EProps::*>::type EPDistanceMap;
typedef boost::property_map<Graph, boost::vertex_index_t>::type VertexIndexMap;



}	// namespace batching_pomp

#endif // BATCHING_POMP_GRAPHTYPES_HPP_