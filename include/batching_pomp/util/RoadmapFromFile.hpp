/***********************************************************************
Copyright (c) 2017, Carnegie Mellon University
All rights reserved.
Authors:  Chris Dellin <cdellin@gmail.com>
          Shushman Choudhury <shushmanchoudhury@gmail.com>

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

#ifndef BATCHING_POMP_UTIL_ROADMAPFROMFILE_HPP_
#define BATCHING_POMP_UTIL_ROADMAPFROMFILE_HPP_

#include <sstream>
#include <iostream>
#include <fstream>
#include <type_traits>
#include <boost/property_map/dynamic_property_map.hpp>
#include <boost/graph/graphml.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <ompl/base/State.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include "rv_state_map_string_adaptor.hpp"

namespace batching_pomp {
namespace util {

//! The map used to decode the .graphml file and populate vertex states.

/// @tparam PropMap The type of property map for vertex states
/// @tparam StateCon The wrapper for the ompl state
template <class PropMap, class StateCon>
class RoadmapFromFilePutStateMap
{
public:
  typedef boost::writable_property_map_tag category;
  typedef typename boost::property_traits<PropMap>::key_type key_type;
  typedef std::string value_type;
  typedef std::string reference;

  const PropMap mPropMap;
  ompl::base::StateSpacePtr mSpace;
  const unsigned int mDim;

  RoadmapFromFilePutStateMap(PropMap _propMap, ompl::base::StateSpacePtr _space, unsigned int _dim)
  : mPropMap{_propMap}
  , mSpace{_space}
  , mDim{_dim}
  {
  }
};

/// Do not allow calling get on this property map
template <class PropMap, class StateCon>
inline std::string
get(const RoadmapFromFilePutStateMap<PropMap,StateCon> & map,
   const typename RoadmapFromFilePutStateMap<PropMap,StateCon>::key_type & k)
{
   abort();
}

/// Convert string representation of vector state space to ompl state
template <class PropMap, class StateCon>
inline void
put(const RoadmapFromFilePutStateMap<PropMap,StateCon> & map,
    const typename RoadmapFromFilePutStateMap<PropMap,StateCon>::key_type & k,
    const std::string representation)
{
  get(map.mPropMap, k).reset(new StateCon(map.mSpace));
  ompl::base::State* ver_state{get(map.mPropMap, k)->state};
  double * values{ver_state->as<ompl::base::RealVectorStateSpace::StateType>()->values};
  std::stringstream ss(representation);
  for (unsigned int ui=0; ui<map.mDim; ui++){
    ss >> values[ui];
  }
}

//! Reads a roadmap encoded as a .graphml file and creates the corresponding Boost Graph

/// Read a graphml file and assign vertices to ompl states
/// Optionally, if file has edges, assign edge distances
/// to be the distance between the states.
/// @tparam Graph The type of boost graph used for the roadmaps
/// @tparam VStateMap The type of boost property map for vertex states
/// @tparam StateCon The wrapper type for an ompl state
/// @tparam EDistance The type of property map for edge lengths
template <class Graph, class VStateMap, class StateCon, class EDistance>
class RoadmapFromFile
{
typedef boost::graph_traits<Graph> GraphTypes;
typedef typename GraphTypes::vertex_descriptor Vertex;
typedef typename GraphTypes::vertex_iterator VertexIter;
typedef typename GraphTypes::edge_descriptor Edge;
typedef typename GraphTypes::edge_iterator EdgeIter;

public:

  const std::string mFilename;

  RoadmapFromFile(
    const ompl::base::StateSpacePtr _space,
    std::string _filename)
  : mSpace(_space)
  , mFilename(_filename)
  , mBounds(0)
  {
    if (mSpace->getType() != ompl::base::STATE_SPACE_REAL_VECTOR) {
      throw std::runtime_error("This only supports real vector state spaces!");
    }

    mDim = mSpace->getDimension();
    mBounds = mSpace->as<ompl::base::RealVectorStateSpace>()->getBounds();
  }

  ~RoadmapFromFile() {}

  void generateVertices(Graph& _roadmap,
                        VStateMap _stateMap)
  {
    std::ifstream fp;
    fp.open(mFilename.c_str());

    boost::dynamic_properties props;
    props.property("state",
      RoadmapFromFilePutStateMap<VStateMap,StateCon>(_stateMap, mSpace, mDim));

    boost::read_graphml(fp, _roadmap, props);
  }

  void generateEdges(Graph& _roadmap,
                     VStateMap _stateMap,
                     EDistance _distanceMap)
  {
    EdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=edges(_roadmap); ei!=ei_end; ++ei)
    {
      ompl::base::State * state1 = get(_stateMap, source(*ei,_roadmap))->state;
      ompl::base::State * state2 = get(_stateMap, target(*ei,_roadmap))->state;
      put(_distanceMap, *ei, mSpace->distance(state1,state2));
    }
  }

private:

  unsigned int mDim;
  ompl::base::RealVectorBounds mBounds;
  const ompl::base::StateSpacePtr mSpace;

};



} // namespace util
} // namespace batching_pomp

#endif //BATCHING_POMP_UTIL_ROADMAPFROMFILE_HPP_