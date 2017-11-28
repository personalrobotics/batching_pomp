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

#ifndef BATCHING_POMP_DUMMY_BATCHING_HPP_
#define BATCHING_POMP_DUMMY_BATCHING_HPP_

#include <ompl/base/StateSpace.h>
#include <ompl/util/Console.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include "batching_pomp/batching/BatchingManager.hpp"

namespace batching_pomp{
namespace batching {

//! Derived class of BatchingManager that implements a single batch search.
  
/// Implements the case where there is a single, reasonably sized roadmap
/// that will be searched by POMP, without any batching. This is
/// to allow a common framework to run POMP on a single roadmap.
class DummyBatching : public BatchingManager
{

public:

  DummyBatching(const ompl::base::StateSpacePtr _space,
                std::string _roadmapFileName
                )
  : BatchingManager(_space,_roadmapFileName)
  {
    BatchingManager::mBatchingType = "dummy";
  }

  //////////////////////////////////////////////////
  /// Overriden methods
  void updateWithNewSolutionCost(double _newSolnCost) override
  {
  }

  void nextBatch(const std::function<bool(const ompl::base::State*)>& _pruneFunction,
                 ompl::NearestNeighbors<Vertex>& _vertexNN) override
  {
    if(BatchingManager::mExhausted){
      OMPL_INFORM("Batching exhausted! No updates with nextBatch!");
      return;
    }

    OMPL_INFORM("Single Batch called!");
    ++BatchingManager::mNumBatches;

    // You know there is only one batch
    // Now remove all invalid vertices
    //BatchingManager::pruneVertices(_pruneFunction,_vertexNN);

    BatchingManager::mExhausted = true;

  }

};



} // namespace batching
} // namespace batching_pomp
#endif // BATCHING_POMP_DUMMY_BATCHING_HPP_
