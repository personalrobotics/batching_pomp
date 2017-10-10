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

#ifndef BATCHING_POMP_UTIL_SELECTOR_HPP_
#define BATCHING_POMP_UTIL_SELECTOR_HPP_

#include <algorithm>
#include <stdexcept>
#include <memory>
#include <functional>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/adjacency_list.hpp>

namespace batching_pomp {
namespace util {

//! A selector assigns an ordering for evaluating edges on a candidate path.

/// Implements various methods of selecting the order of edges
/// to check for lazy evaluation of a candidate path.
/// @tparam Graph The type of boost graph used for the roadmaps
template<class Graph>
class Selector
{

typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::edge_descriptor Edge;

public:

  /// \param[in] _type The kind of edge selector to use.
  ///             Options are normal, alternate, failfast and maxinf
  Selector(const std::string& _type)
  : mType{_type}
  {
    if(mType == "normal") {
      
      mSelectEdgeFnPtr = std::bind(&Selector<Graph>::selectEdgesNormal,this,std::placeholders::_1,std::placeholders::_2);
    }
    else if(mType == "alternate") {
      mSelectEdgeFnPtr = std::bind(&Selector<Graph>::selectEdgesAlternate,this,std::placeholders::_1,std::placeholders::_2);
    }
    else if(mType == "failfast") {
      mSelectEdgeFnPtr = std::bind(&Selector<Graph>::selectEdgesFailFast,this,std::placeholders::_1,std::placeholders::_2);
    }
    else if(mType == "maxinf") {
      mSelectEdgeFnPtr = std::bind(&Selector<Graph>::selectEdgesMaxInf,this,std::placeholders::_1,std::placeholders::_2);
    }
    else {
      throw std::invalid_argument("Invalid selector type specified - "+mType+"!");
    }
  }

  std::string getType()
  {
    return mType;
  }

  /// Call edge selector function pointer
  std::vector<Edge> selectEdges(const Graph& g, const std::vector<Edge>& epath) const
  {
    return mSelectEdgeFnPtr(g,epath);
  }


private:
  // For a pair of edge and probability of collision
  using EdgePair = std::pair<double,Edge>;

  // Just return in order of edges on path
  std::vector<Edge> selectEdgesNormal(const Graph& g, const std::vector<Edge>& epath) const
  {
    return epath;
  }

  // Return reverse order of edges
  std::vector<Edge> selectEdgesAlternate(const Graph& g, const std::vector<Edge>& epath) const
  {
    std::vector<Edge> reverseEpath(epath);
    std::reverse(reverseEpath.begin(),reverseEpath.end());
    return reverseEpath;
  }

  // Return in order of highest probability of collision
  std::vector<Edge> selectEdgesFailFast(const Graph& g, const std::vector<Edge>& epath) const
  {
    std::vector < EdgePair > edgeMeasList;
    edgeMeasList.reserve(epath.size());
    for(Edge e : epath) {
      edgeMeasList.push_back(std::make_pair(g[e].collMeasure,e));
    }

    std::sort(edgeMeasList.rbegin(), edgeMeasList.rend());

    std::vector<Edge> orderedEdges;
    orderedEdges.reserve(edgeMeasList.size());
    for(EdgePair ep : edgeMeasList) {
      orderedEdges.push_back(ep.second);
    }

    return orderedEdges;
  }

  // Return in order or most uncertain (closest to 0.5 probability)
  std::vector<Edge> selectEdgesMaxInf(const Graph& g, const std::vector<Edge>& epath) const
  {
    std::vector < EdgePair > edgeProbList;
    edgeProbList.reserve(epath.size());
    for(Edge e : epath) {
      edgeProbList.push_back( std::make_pair(std::fabs(std::exp(-g[e].collMeasure)-0.5),e) );
    }

    std::sort(edgeProbList.begin(), edgeProbList.end());

    std::vector<Edge> orderedEdges;
    orderedEdges.reserve(edgeProbList.size());
    for(EdgePair ep : edgeProbList) {
      orderedEdges.push_back(ep.second);
    }

    return orderedEdges;
  }

  std::string mType;
  std::function< std::vector<Edge>(const Graph& g, const std::vector<Edge>& epath) > mSelectEdgeFnPtr;

};

} // namespace util
} // namespace batching_pomp

#endif // BATCHING_POMP_UTIL_SELECTOR_HPP_