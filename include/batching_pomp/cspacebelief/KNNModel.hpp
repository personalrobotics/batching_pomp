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

#ifndef BATCHING_POMP_CSPACEBELIEF_KNNMODEL_HPP_
#define BATCHING_POMP_CSPACEBELIEF_KNNMODEL_HPP_

#include <functional>
#include <exception>
#include <memory>
#include <ompl/datastructures/NearestNeighbors.h>
#include <ompl/datastructures/NearestNeighborsLinear.h>
#include <ompl/datastructures/NearestNeighborsGNAT.h>
#include <Eigen/Dense>
#include "Model.hpp"
#include "BeliefPoint.hpp"

namespace batching_pomp {
namespace cspacebelief {

/// Implements a C-Space Belief Model using k-Nearest Neighbour lookup
/// And weighted averaging over neighbouring belief points
/// The weights are inversely proportional to the distance

class KNNModel : public virtual Model<BeliefPoint>
{
public:

  /// \param[in] _KNN The value of k for KNN lookup
  /// \param[in] _supportThreshold The distance threshold to consider nearest neighbors for
  /// \param[in] _distanceFunction The ompl NearestNeighbours Distance Function to set 
  KNNModel(size_t _KNN, double _supportThreshold,
           double _prior, double _priorWeight,
           const std::function<double(const BeliefPoint&, const BeliefPoint&)>& _distanceFunction)
  : mNumPoints{0}
  , mKNN{_KNN}
  , mSupportThreshold{_supportThreshold}
  , mPrior{_prior}
  , mPriorWeight{_priorWeight}
  , mDistanceFunction{_distanceFunction}
  {
    mBeliefPointNN.reset(new ompl::NearestNeighborsGNAT<BeliefPoint>());
    mBeliefPointNN->setDistanceFunction(mDistanceFunction);
  }

  //////////////////////////////////////////////////
  // Setters
  void setPrior(double _prior)
  {
    mPrior = _prior;
  }

  void setPriorWeight(double _priorWeight)
  {
    mPriorWeight = _priorWeight;
  }

  void setKNN(size_t _knn)
  {
    mKNN = _knn;
  }

  void setSupportThreshold(double _supportThreshold)
  {
    mSupportThreshold = _supportThreshold;
  }

  //////////////////////////////////////////////////
  // Overriden methods
  void addPoint(const BeliefPoint& data) override
  {
    mBeliefPointNN->add(data);
    mNumPoints++;
  }

  void removePoint(const BeliefPoint& data) override
  {
    mBeliefPointNN->remove(data);
    mNumPoints--;
  }

  double estimate(const BeliefPoint& query) const override
  {
    if(mNumPoints == 0) {
      return mPrior;
    }

    // If fewer than k thus far, take all points
    size_t knn = std::min(mKNN , mNumPoints);

    std::vector<BeliefPoint> neighboursVect;
    mBeliefPointNN->nearestK(query,knn,neighboursVect);

    if(knn != neighboursVect.size()) {
      throw std::runtime_error("Model could not return the expected number of neighbours");
    }

    Eigen::VectorXd weights(knn);
    Eigen::VectorXd values(knn);

    for(size_t i = 0; i < knn; i++) {
      double distance = mDistanceFunction(query, neighboursVect[i]);

      if(distance < mSupportThreshold) {
        weights[i] = 1.0/distance;
        values[i] = neighboursVect[i].getValue(); // Second element of pair is value
      }
    }

    double result = weights.dot(values) / weights.sum();

    // Do (result + pw*p) / (1 + pw) for smoothing by prior
    result = (result + mPriorWeight*mPrior) / (1 + mPriorWeight);

    return result;
  }

private:
    
  size_t mNumPoints;
  size_t mKNN;
  double mSupportThreshold;
  double mPrior;
  double mPriorWeight;
  std::unique_ptr<ompl::NearestNeighborsGNAT<BeliefPoint>> mBeliefPointNN;
  std::function<double(const BeliefPoint&, const BeliefPoint&)> mDistanceFunction;
};

} //namespace cspacebelief
} //namespace batching_pomp

#endif //BATCHING_POMP_CSPACEBELIEF_KNNMODEL_HPP_