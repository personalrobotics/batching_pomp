/***********************************************************************
Copyright (c) 2016, Carnegie Mellon University
All rights reserved.
Authors: Christopher Dellin <cdellin@gmail.com>
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
#ifndef BATCHING_POMP_UTIL_RVSTATEMAP_H
#define BATCHING_POMP_UTIL_RVSTATEMAP_H

#include <boost/property_map/dynamic_property_map.hpp>
#include <ompl/base/State.h>
#include <ompl/base/StateSpace.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>


namespace batching_pomp {
namespace util {

inline void stringify_from_x(std::string & repr, const double & val)
{
    char buf[2048];
    // invariant: min DOESNT WORK, max DOES WORK
    // validate invariants
    int min = 0;
    sprintf(buf, "%.*f", min, val);
    if (val == strtod(buf,0))
    {
        repr = std::string(buf);
        return;
    }
    // what is it at 1?
    sprintf(buf, "%.*f", 1, val);
    int max = sizeof(buf)-strlen(buf);
    sprintf(buf, "%.*f", max, val);
    if (val != strtod(buf,0))
    {
        printf("stringify_from_x invariant failed!\n");
        abort();
    }
    // binary search
    for (;;)
    {
        int diff = max - min;
        if (diff == 1)
            break;
        int test = min + diff/2;
        sprintf(buf, "%.*f", test, val);
        if (val == strtod(buf,0))
            max = test;
        else
            min = test;
    }
    sprintf(buf, "%.*f", max, val);
    repr = std::string(buf);
    return;
}

/*! \brief Treats OMPL state space values as space-delimited strings
* Allocates a new state while putting, but does not free existing ones.
*/

template <class StateMap>
class rvstate_map_string_adaptor
{
public:
    typedef boost::read_write_property_map_tag category;
    typedef typename boost::property_traits<StateMap>::key_type key_type;
    typedef std::string value_type;
    typedef std::string reference;
    const StateMap state_map;
    ompl::base::RealVectorStateSpace * rvspace;
    const unsigned int dim;
    rvstate_map_string_adaptor(StateMap state_map, ompl::base::RealVectorStateSpace * rvspace):
      state_map(state_map), rvspace(rvspace), dim(rvspace->getDimension())
    {
    }
};

template <class StateMap>
inline std::string
get(const rvstate_map_string_adaptor<StateMap> & adaptor,
   const typename rvstate_map_string_adaptor<StateMap>::key_type & k)
{
    ompl::base::RealVectorStateSpace::StateType * rvstate
        = (ompl::base::RealVectorStateSpace::StateType *)get(adaptor.state_map, k);
    if (!rvstate)
        return std::string();
    std::string s;
    for (unsigned int ui=0; ui<adaptor.dim; ui++)
    {
        if (ui) s += " ";
        std::string component_repr;
        util::stringify_from_x(component_repr, rvstate->values[ui]);
        s += component_repr;
    }
    return s;
}

template <class StateMap>
inline void
put(const rvstate_map_string_adaptor<StateMap> & adaptor,
   const typename rvstate_map_string_adaptor<StateMap>::key_type & k,
   const std::string s)
{
   ompl::base::RealVectorStateSpace::StateType * rvstate;
   if (s.length() == 0)
   {
      rvstate = 0;
   }
   else
   {
      rvstate = (ompl::base::RealVectorStateSpace::StateType *)adaptor.rvspace->allocState();
      std::stringstream ss(s);
      for (unsigned int ui=0; ui<adaptor.dim; ui++)
         ss >> rvstate->values[ui];
   }
   put(adaptor.state_map, k, rvstate);
}

template <class StateMap>
rvstate_map_string_adaptor<StateMap>
make_rvstate_map_string_adaptor(StateMap state_map, ompl::base::RealVectorStateSpace * rvspace)
{
   return rvstate_map_string_adaptor<StateMap>(state_map, rvspace);
}


}
}

#endif