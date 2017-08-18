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

#ifndef BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_
#define BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_

#include <ompl/base/State.h>
#include <Eigen/Dense>
#include <functional>

namespace batching_pomp {
namespace cspacebelief {

//! An example datatype of the configuration space belief model.

/// Consists of the state and the result
/// of the collision check of that state. The assumption is
/// that every belief model will at a minimum require the value
/// of the checked configuration and its result.
struct BeliefPoint {

  Eigen::VectorXd stateValues;
  double value;

  BeliefPoint():value(0.0){}
  BeliefPoint(const ompl::base::State* _state, 
              unsigned int _dims, double _val)
  :value{_val}
  {
    double *values = _state->as<ompl::base::RealVectorStateSpace::StateType>()->values;
    stateValues = Eigen::Map<Eigen::VectorXd>(values,_dims);
  }

  double getValue()
  {
    return value;
  }

  bool operator==(const BeliefPoint& bp) const
  {
    return (stateValues==bp.stateValues);
  }

  bool operator!=(const BeliefPoint& bp) const
  {
    return (stateValues!=bp.stateValues);
  }
};

} //namespace cspacebelief
} //namespace batching_pomp

#endif //BATCHING_POMP_CSPACEBELIEF_BELIEFPOINT_HPP_