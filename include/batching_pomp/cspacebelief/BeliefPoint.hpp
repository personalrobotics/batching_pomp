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

#ifndef BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_
#define BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_

#include <ompl/base/State.h>
#include <functional>

namespace batching_pomp {
namespace cspacebelief {

struct BeliefPoint {

  const ompl::base::State* state;
  double value;

  BeliefPoint():value(0.0){}
  BeliefPoint(ompl::base::State* _state, double _val)
    : state(_state)
    , value(_val) {}

  double getValue()
  {
    return value;
  }

  // TODO : Are these operator overloads correct?
  bool operator==(const BeliefPoint& bp) const
  {
    return (state==bp.state);
  }

  bool operator!=(const BeliefPoint& bp) const
  {
    return (state!=bp.state);
  }
};

/// Computes (symmetric) distance metric between two belief-point instances
/// \param[in] _distFun The distance function between underlying ompl state instances
/// \param[in] _bp1, _bp2 The two belief point instances
/// \return The distance between two belief points
double beliefDistanceFunction(
  std::function<double(const ompl::base::State*, const ompl::base::State*)>& _distFun,
  const BeliefPoint& _bp1, const BeliefPoint& _bp2)
{
  return _distFun(_bp1.state, _bp2.state);
}

} //namespace cspacebelief
} //namespace batching_pomp

#endif //BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_