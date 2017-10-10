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
#ifndef BATCHING_POMP_BELIEF_RESAMPLING_BELIEFINFORMEDRESAMPLER_HPP_
#define BATCHING_POMP_BELIEF_RESAMPLING_BELIEFINFORMEDRESAMPLER_HPP_

#include <functional>

#include <ompl/base/StateSpace.h>
#include <ompl/datastructures/NearestNeighbors.h>

#include "batching_pomp/cspacebelief.hpp"
#include "batching_pomp/GraphTypes.hpp"
#include "batching_pomp/BatchingPOMP.hpp"

namespace batching_pomp {
namespace belief_resampling {

class BeliefInformedResampler
{
public:

BeliefInformedResampler(BatchingPOMP& _planner)
: mPlanner(_planner)
, mBeliefModel{_planner.mBeliefModel}
, mSpace{_planner.mSpace}
, mCurrentRoadmap(*(_planner.g))
{
}

virtual ~BeliefInformedResampler() = default;

void setNumberOfSamples(unsigned int _nSamples)
{
	mBatchParams.first = _nSamples;
}

void setConnectionRadius(double _radius)
{
	mBatchParams.second = _radius;
}

void setBatchParams(unsigned int _nSamples, double _radius)
{
	setNumberOfSamples(_nSamples);
	setConnectionRadius(_radius);
}

virtual void updateRoadmap() = 0;


protected:

	/// The planner object itself
	BatchingPOMP& mPlanner;

	/// The pointer to the C-space belief model instance to be used by the sampler
	std::shared_ptr< cspacebelief::Model<cspacebelief::BeliefPoint> > mBeliefModel;

	/// The pointer to the OMPL state space 
	const ompl::base::StateSpacePtr mSpace;

	/// The underlying roadmap of the planner
	Graph& mCurrentRoadmap;

	/// Pair of samples and radius
	std::pair<unsigned int, double> mBatchParams;

};

} //namespace cspacebelief
} //namespace batching_pomp

#endif