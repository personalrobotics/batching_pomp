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

#include "batching_pomp/include/BatchingPOMP.hpp"

namespace batching_pomp{


/// Default constructor for CPP-level usage
BatchingPOMP::BatchingPOMP(const ompl::base::SpaceInformationPtr & si,
                           std::unique_ptr<BatchingManager<Graph,VPStateMap,StateCon,EPDistanceMap>> _batchingPtr,
                           std::unique_ptr< Model<BeliefPoint> _beliefModel,
                           std::unique_ptr<Selector<Graph>> _selector,
                           double _searchInflationFactor,
                           double _decrement,
                           double _startGoalRadius,
                           double _checkRadius)
: ompl::base::Planner(si,"BatchingPOMP")
, mSpace(si->getStateSpace()),
, mBatchingPtr{std::move(_batchingPtr)}
, mBeliefModel{std::move(_beliefModel)}
, mSelector{std::move(_selector)}
, mSearchInflationFactor{_searchInflationFactor}
, mDecrement{_decrement}
, mStartGoalRadius{_startGoalRadius}
, mIsInitSearchBatch{false}
, mIsPathfound{false}
, mCheckRadius{_checkRadius}
, mBestCost{std::numeric_limits<double>::max()}
, mNumEdgeChecks{0u}
, mNumCollChecks{0u}
, mNumSearches{0u}
, mLookupTime{0.0}
{
  
}





} // namespace batching_pomp