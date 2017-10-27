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
#ifndef BATCHING_POMP_BELIEF_RESAMPLING_SGDBASEDRESAMPLER_HPP_
#define BATCHING_POMP_BELIEF_RESAMPLING_SGDBASEDRESAMPLER_HPP_

#include <cmath>
#include <exception>
#include <random>
#include <thread>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <ompl/base/StateSpace.h>
#include <ompl/base/ScopedState.h>
#include <ompl/base/spaces/RealVectorStateSpace.h>
#include <ompl/util/Exception.h>
#include <Eigen/Dense>

#include "batching_pomp/GraphTypes.hpp"
#include "batching_pomp/cspacebelief.hpp"
#include "BeliefInformedResampler.hpp"
#include "batching_pomp/BatchingPOMP.hpp"

namespace batching_pomp {
namespace belief_resampling {


class SGDBasedResampler : public virtual BeliefInformedResampler
{

public:

  SGDBasedResampler(BatchingPOMP& _planner,
                    unsigned int _numPerturbations,
                    double _dAlpha,
                    double _perturbSize = 0.1,
                    double _probThreshold = 1.0)
  : BeliefInformedResampler(_planner)
  , mFullRoadmap(*(_planner.full_g))
  , mCurrVertexNN(*(_planner.mVertexNN))
  , mStartVertex{_planner.getStartVertex()}
  , mGoalVertex{_planner.getGoalVertex()}
  , mNumPerturbations{_numPerturbations}
  , mPerturbSize{_perturbSize}
  , mProbThreshold{_probThreshold}
  , mSimpleFlag{false}
  , mSuccPerturbations{0u}
  {
    unsigned int nAlphas = static_cast<unsigned int>(1.0/_dAlpha) + 1u;
    alphaMinVector.resize(nAlphas);
    alphaMinIdxVector.resize(nAlphas);
    alphaSecondMinVector.resize(nAlphas);
    boost::tie(mCurrVertex,mLastVertex) = vertices(mFullRoadmap);
    std::random_device rd;
    mGen.seed(rd());
  }

  void setSimpleFlag()
  {
    mSimpleFlag = true;
  }

  /// Add the next n samples from list subject to thresholding conditions
  void addInitialBatch()
  {
    // Iterate and add to current roadmap
    unsigned int numVerticesAdded{0u};
    std::vector<Vertex> vertex_vector(mBatchParams.first);

    while(numVerticesAdded < mBatchParams.first)
    {
      cspacebelief::BeliefPoint query(mFullRoadmap[*mCurrVertex].v_state->state, 
                        mSpace->getDimension(), -1.0);

      // TODO - Currently coll check thresholding here; update
      if(mPlanner.getSpaceInformation()->isValid(mFullRoadmap[*mCurrVertex].v_state->state) &&
         mBeliefModel->estimate(query) < mProbThreshold) {
        Vertex newVertex{boost::add_vertex(mCurrentRoadmap)};

        mCurrentRoadmap[newVertex].v_state = mFullRoadmap[*mCurrVertex].v_state;
        vertex_vector[numVerticesAdded++] = newVertex;
      }
      ++mCurrVertex;
    }

    if(numVerticesAdded > 0u) {
      vertex_vector.resize(numVerticesAdded);
      mCurrVertexNN.add(vertex_vector);
    }
    else{
      throw ompl::Exception("No vertices satisfied probability threshold!");
    }

    mNumCurrVertices = boost::num_vertices(mCurrentRoadmap);

    // Define edges based on mutual distance and compute edge weights
    VertexIter vi, vi_end;
    for(boost::tie(vi,vi_end)=vertices(mCurrentRoadmap); vi != vi_end; ++vi) {
      std::vector<Vertex> vertexNbrs;
      mCurrVertexNN.nearestR(*vi,mBatchParams.second,vertexNbrs);

      for(Vertex nbr : vertexNbrs) {
        if(nbr == *vi)
          continue;
        if(!boost::edge(*vi, nbr, mCurrentRoadmap).second) {
          std::pair<Edge,bool> new_edge = boost::add_edge(*vi,nbr,mCurrentRoadmap);
          mCurrentRoadmap[new_edge.first].distance = mPlanner.vertexDistFun(*vi,nbr);
          mCurrentRoadmap[new_edge.first].blockedStatus = UNKNOWN;
          mCurrentRoadmap[new_edge.first].hasPoints = false;
          mCurrentRoadmap[new_edge.first].collMeasure = mPlanner.computeEdgeCollisionMeasureNoStates(*vi, nbr);
        }
      }
    }
  }


  void computeInitialScore()
  {
    vertexImportance.resize(mNumCurrVertices,0.0);

    mInfiniteCostEdges.clear();

    // RUN 4 DIJKSTRA SEARCHES WITH (Start,Goal) X (M,L)

    // First for M_to_come and L_to_come from start
    boost::dijkstra_shortest_paths(mCurrentRoadmap, mStartVertex,
                                   boost::make_assoc_property_map(M_preds),
                                   boost::make_assoc_property_map(M_to_come),
                                   boost::get(&EProps::collMeasure,mCurrentRoadmap),
                                   boost::get(boost::vertex_index, mCurrentRoadmap),
                                   std::less<double>(),
                                   boost::closed_plus<double>(std::numeric_limits<double>::max()), 
                                   std::numeric_limits<double>::max(),
                                   double(),
                                   boost::make_dijkstra_visitor(boost::null_visitor()));

    boost::dijkstra_shortest_paths(mCurrentRoadmap, mStartVertex,
                                   boost::make_assoc_property_map(L_preds),
                                   boost::make_assoc_property_map(L_to_come),
                                   boost::get(&EProps::distance,mCurrentRoadmap),
                                   boost::get(boost::vertex_index, mCurrentRoadmap),
                                   std::less<double>(),
                                   boost::closed_plus<double>(std::numeric_limits<double>::max()), 
                                   std::numeric_limits<double>::max(),
                                   double(),
                                   boost::make_dijkstra_visitor(boost::null_visitor()));

    boost::dijkstra_shortest_paths(mCurrentRoadmap, mGoalVertex,
                                   boost::make_assoc_property_map(M_succs),
                                   boost::make_assoc_property_map(M_to_go),
                                   boost::get(&EProps::collMeasure,mCurrentRoadmap),
                                   boost::get(boost::vertex_index, mCurrentRoadmap),
                                   std::less<double>(),
                                   boost::closed_plus<double>(std::numeric_limits<double>::max()), 
                                   std::numeric_limits<double>::max(),
                                   double(),
                                   boost::make_dijkstra_visitor(boost::null_visitor()));

    boost::dijkstra_shortest_paths(mCurrentRoadmap, mGoalVertex,
                                   boost::make_assoc_property_map(L_succs),
                                   boost::make_assoc_property_map(L_to_go),
                                   boost::get(&EProps::distance,mCurrentRoadmap),
                                   boost::get(boost::vertex_index, mCurrentRoadmap),
                                   std::less<double>(),
                                   boost::closed_plus<double>(std::numeric_limits<double>::max()), 
                                   std::numeric_limits<double>::max(),
                                   double(),
                                   boost::make_dijkstra_visitor(boost::null_visitor()));


    // PRINT MAPS AND ENSURE ALL THERE

    // FOR EACH ALPHA, GET VECTORS OVER VERTICES AND COMPUTE SCORES
    // ALSO COMPUTE ALPHA SECOND MIN VECTOR
    unsigned int nAlphas = static_cast<unsigned int>(alphaMinVector.size());

    for(unsigned int i = 0u; i < nAlphas; i++) {

      double currAlpha{static_cast<double>(i)/(nAlphas-1)};

      std::vector<double> currAlphaScores(mNumCurrVertices);

      for (auto verValPair : L_to_come)
      {
        currAlphaScores[verValPair.first] = currAlpha*(verValPair.second + L_to_go[verValPair.second])
          + (1.0 - currAlpha)*(M_to_come[verValPair.first] + M_to_go[verValPair.first]);
      }

      // Use custom min and second min to populate alpha vectors
      double min = std::numeric_limits<double>::infinity();
      unsigned int min_idx{0u};
      double secondMin = std::numeric_limits<double>::infinity();

      for(unsigned int j=0u; j< mNumCurrVertices; j++) {
        if(currAlphaScores[j] < min) {
          secondMin = min;
          min = currAlphaScores[j];
          min_idx = j;
        }
        else if(currAlphaScores[j] < secondMin && currAlphaScores[j] > min) {
          secondMin = currAlphaScores[j];
        }
      }

      alphaMinVector[i] = min;
      alphaMinIdxVector[i] = min_idx;
      alphaSecondMinVector[i] = secondMin;
    
      // Now iterate over vertices and assign importance
      for(unsigned int j=0u; j< mNumCurrVertices; j++) {
        if(std::fabs(currAlphaScores[j] - min) < std::numeric_limits<double>::epsilon()) {
          vertexImportance[j] += secondMin - min;
        }
      }
    }

    mCurrRoadmapScore = getScoreFromAlphaVector(alphaMinVector,alphaMinIdxVector);
    std::cout<<"currentRoadmapScore - "<<mCurrRoadmapScore<<std::endl;
  }

  /// The score is basically the sum of the
  /// vector divided by the number of unique values (paths)
  double getScoreFromAlphaVector(const std::vector<double>& _alphaMinVector,
                                 const std::vector<unsigned int>& _alphaMinIdxVector) const
  {
    unsigned int currIdx{_alphaMinIdxVector[0]};
    double score{_alphaMinVector[0]};
    unsigned int numUniqueVals{1u};
    unsigned int nAlphas = static_cast<unsigned int>(alphaMinVector.size());

    for(unsigned int i = 1u; i < nAlphas; i++) {
      score += _alphaMinVector[i];
      if(_alphaMinIdxVector[i] != currIdx) {
        numUniqueVals++;
        currIdx = _alphaMinIdxVector[i];
      }
    }

    score /= numUniqueVals;
    return score;
  }


  bool perturbVertexNaive(const Vertex& u, StateConPtr& perturbedState) const
  {
    //StateConPtr perturbedState(std::make_shared<StateCon>(mSpace));

    ompl::base::RealVectorStateSampler rvSampler(mSpace.get());

    rvSampler.sampleUniformNear(perturbedState->state, mCurrentRoadmap[u].v_state->state, mPerturbSize);

    if(mPlanner.getSpaceInformation()->isValid(perturbedState->state)) {
      return true;
    }

    return false;
  }



  bool perturbVertexApproxGrad(const Vertex& u, StateConPtr& perturbedState) const
  {
    // Iterate through neighbors of u and compute L and M and sum
    double vert_LM_pdt{0.0};
    unsigned int idx{0u};
    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei) {
      vert_LM_pdt += mCurrentRoadmap[*ei].distance * mCurrentRoadmap[*ei].collMeasure;
      idx++;
    }

    // Get estimate of vertex
    unsigned int nDims{mSpace->getDimension()};
    cspacebelief::BeliefPoint query(mCurrentRoadmap[u].v_state->state, 
                                    nDims, -1.0);
    double vertLogProb{-std::log(mBeliefModel->estimate(query))};

    double GRAD_EPSILON{0.001};

    Eigen::VectorXd grads(nDims);
    //std::cout<<"Coords - ";
    for(unsigned int d=0u; d < nDims; d++)
    {
      StateConPtr oneDimPerturbed(std::make_shared<StateCon>(mSpace));
      mSpace->copyState(oneDimPerturbed->state, mCurrentRoadmap[u].v_state->state);

      // Perturb along dimension
      double* values = oneDimPerturbed->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      //std::cout<<values[d]<<" ";
      values[d] = values[d] + GRAD_EPSILON;

      // Now compute new edge lengths and measures and compute difference
      double new_LM_pdt{0.0};

      for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei) {
        Vertex v{boost::target(*ei,mCurrentRoadmap)};

        // compute perturbed length and measure
        double new_length{mSpace->distance(oneDimPerturbed->state, mCurrentRoadmap[v].v_state->state)};

        cspacebelief::BeliefPoint per_query(oneDimPerturbed->state, 
                        mSpace->getDimension(), -1.0);
        double newLogProb{-std::log(mBeliefModel->estimate(per_query))};
        double new_meas = mCurrentRoadmap[*ei].collMeasure + newLogProb - vertLogProb;
        new_LM_pdt += new_length*new_meas;
      }

      grads[d] = (new_LM_pdt - vert_LM_pdt) / GRAD_EPSILON;
    }

    // Normalize vector
    grads.stableNormalize();
    //std::cout<<std::endl<<"GRADS - "<<grads<<std::endl;

    if(isnan(grads[0]) || isnan(grads[1])){
      return false;
    }

    // Now what????
    // Do backtracking line search till collision-free?
    double stepSize{mPerturbSize};
    StateConPtr perturbedStateCopy(std::make_shared<StateCon>(mSpace));
    while(stepSize > 0.01)
    {
      mSpace->copyState(perturbedStateCopy->state, mCurrentRoadmap[u].v_state->state);
      double* values = perturbedStateCopy->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      for(unsigned int d=0u; d < nDims; d++)
      {
        values[d] = values[d] + stepSize*grads[d];
      }

      if( mSpace->satisfiesBounds(perturbedStateCopy->state) &&
            mPlanner.getSpaceInformation()->isValid(perturbedStateCopy->state)) {
        break;
      }

      stepSize /= 2.0;
    }

    //std::cout<<stepSize<<std::endl;

    if(stepSize > 0.01) {
      mSpace->copyState(perturbedState->state,perturbedStateCopy->state);
      return true;
    }

    return false;

  }


  double implementPerturbation(Vertex& u, const StateConPtr& _perturbedState)
  {
    cspacebelief::BeliefPoint query(_perturbedState->state, 
                        mSpace->getDimension(), -1.0);
    if(mBeliefModel->estimate(query) > mProbThreshold) {
      return 0.0;
    }
    // Use this as a proxy for new vertex
    Vertex tempNewVertex{boost::add_vertex(mCurrentRoadmap)};
    
    //mCurrentRoadmap[tempNewVertex].v_state = _perturbedState;
    mCurrentRoadmap[tempNewVertex].v_state = std::make_shared<StateCon>(mSpace);
    mSpace->copyState(mCurrentRoadmap[tempNewVertex].v_state->state, _perturbedState->state);

    // To recompute old info
    std::map<Vertex,std::pair<double,double>> newNbrVertexWeights;

    // First iterate through current edges and make inf cost those that are invalidated
    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei) {

      double edgeDist{mPlanner.vertexDistFun(boost::target(*ei,mCurrentRoadmap), tempNewVertex)};

      

      // Edge invalid in temp new graph
      if(edgeDist > mBatchParams.second) { 
        newNbrVertexWeights.insert(std::make_pair(boost::target(*ei,mCurrentRoadmap),
          std::make_pair(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity())));
      }
      else {
        // Store modified weights
        newNbrVertexWeights.insert(std::make_pair(boost::target(*ei,mCurrentRoadmap),
          std::make_pair(edgeDist, mPlanner.computeEdgeCollisionMeasureNoStates(
            boost::target(*ei,mCurrentRoadmap) , tempNewVertex))));
      }
    }

    std::vector<Vertex> tempVertexNbrs;
    mCurrVertexNN.nearestR(tempNewVertex,mBatchParams.second,tempVertexNbrs);

    double min_M_tc{std::numeric_limits<double>::infinity()};
    double min_M_tg{std::numeric_limits<double>::infinity()};
    double min_L_tc{std::numeric_limits<double>::infinity()};
    double min_L_tg{std::numeric_limits<double>::infinity()};
    Vertex M_pred,M_succ,L_pred,L_succ;

    for(Vertex tvnbr : tempVertexNbrs) {

      std::pair<double,double> weights{std::make_pair(mPlanner.vertexDistFun(tempNewVertex,tvnbr), 
            mPlanner.computeEdgeCollisionMeasureNoStates(tempNewVertex,tvnbr))};

      if(tvnbr != u && !boost::edge(u,tvnbr,mCurrentRoadmap).second) {
        newNbrVertexWeights.insert(std::make_pair(tvnbr,weights));
      }

      //Either way, update the M_tc,M_tg  OF NEW VERTEX etc
      if(L_to_come[tvnbr] + weights.first < min_L_tc) {
        min_L_tc = L_to_come[tvnbr] + weights.first;
        L_pred = tvnbr;
      }
      if(L_to_go[tvnbr] + weights.first < min_L_tg) {
        min_L_tg = L_to_go[tvnbr] + weights.first;
        L_succ = tvnbr;
      }
      if(M_to_come[tvnbr] + weights.second < min_M_tc) {
        min_M_tc = M_to_come[tvnbr] + weights.second;
        M_pred = tvnbr;
      }
      if(M_to_go[tvnbr] + weights.second < min_M_tg) {
        min_M_tg = M_to_go[tvnbr] + weights.second;
        M_succ = tvnbr;
      }
    }

    typedef std::tuple<double,Vertex,double,Vertex,double,Vertex,double,Vertex> vertexInfoTuple;

    //TODO : Just store the new indices for alphaMin etc. and then update in place at end?

    std::vector<double> tempAlphaVector(alphaMinVector);
    std::vector<unsigned int> tempAlphaIdxVector(alphaMinIdxVector);
    std::vector<double> tempSecondAlphaVector(alphaSecondMinVector);
    std::map<Vertex, double> tempVertexImportanceMap;
    std::map<Vertex, vertexInfoTuple> tempVertexInfoMap;

    // Only vertices below can affect the score
    for(auto nvNbrWts : newNbrVertexWeights) {
      // For each tvnbr, recompute best M_tg, M_tc etc.
      // IF either is better than their current, put in map
      Vertex nvnbr{nvNbrWts.first};
      std::pair<double,double> weights{nvNbrWts.second};

      double v_Min_L_tc{std::numeric_limits<double>::infinity()};
      double v_Min_M_tc{std::numeric_limits<double>::infinity()};
      double v_Min_L_tg{std::numeric_limits<double>::infinity()};
      double v_Min_M_tg{std::numeric_limits<double>::infinity()};
      Vertex v_L_pred, v_M_pred, v_L_succ, v_M_succ;

      // If its predecessor OR successor is currently u, reselect best
      // ELSE just compare with current
      if(L_preds[nvnbr] == u || M_preds[nvnbr] == u
         || L_succs[nvnbr] == u || L_preds[nvnbr]==u) {

        OutEdgeIter ei, ei_end;

        for(boost::tie(ei,ei_end)=out_edges(nvnbr,mCurrentRoadmap); ei != ei_end; ++ei) {

          Vertex candidate{boost::target(*ei,mCurrentRoadmap)};
          double length, meas;
          double ltc,mtc,ltg,mtg;

          if(candidate == u) {
            // Re-assign weights to temporary
            length = weights.first;
            meas = weights.second;
            ltc = min_L_tc;
            mtc = min_M_tc;
            ltg = min_L_tg;
            mtg = min_M_tg;
          }
          else{
            length = mCurrentRoadmap[*ei].distance;
            meas = mCurrentRoadmap[*ei].collMeasure;
            ltc = L_to_come[candidate];
            mtc = M_to_come[candidate];
            ltg = L_to_go[candidate];
            mtg = M_to_go[candidate];
          }

          // Update to-go/to-come temporarily
          if(ltc + length < v_Min_L_tc){
            v_Min_L_tc = ltc + length;
            v_L_pred = candidate;
          }
          if(mtc + meas < v_Min_M_tc){
            v_Min_M_tc = mtc + meas;
            v_M_pred = candidate;
          }
          if(ltg + length < v_Min_L_tg){
            v_Min_L_tg = ltg + length;
            v_L_succ = candidate;
          }
          if(mtg + meas < v_Min_M_tg){
            v_Min_M_tg = mtg + meas;
            v_M_succ = candidate;
          }
        }
      }
      else{
        //Just update to tempNewVertex if it can be reached more easily
        if(min_L_tc + weights.first < L_to_come[nvnbr]) {
          v_Min_L_tc = min_L_tc + weights.first;
          v_L_pred = u;
        }
        if(min_M_tc + weights.second < M_to_come[nvnbr]) {
          v_Min_M_tc = min_M_tc + weights.second;
          v_M_pred = u;
        }
        if(min_L_tg + weights.first < L_to_go[nvnbr]) {
          v_Min_L_tg = min_L_tg + weights.first;
          v_L_succ = u;
        }
        if(min_M_tg + weights.second < M_to_go[nvnbr]) {
          v_Min_M_tg = min_M_tg + weights.second;
          v_M_succ = u;
        }
      }


      // Now test if nvnbr's M/L has changed
      // And if so, put in map and compute effect on score
      if(v_Min_L_tc + v_Min_L_tg < L_to_come[nvnbr] + L_to_go[nvnbr]
        || v_Min_M_tc + v_Min_M_tg < M_to_come[nvnbr] + M_to_go[nvnbr])
      {

        vertexInfoTuple nvnbr_tuple = 
          std::make_tuple(v_Min_L_tc, v_L_pred, v_Min_M_tc, v_M_pred,
                          v_Min_L_tg, v_L_succ, v_Min_M_tg, v_M_succ);
        tempVertexInfoMap.insert(std::make_pair(nvnbr, nvnbr_tuple));
        tempVertexImportanceMap.insert(std::make_pair(nvnbr, vertexImportance[nvnbr]));
      }
    }


    // Now run through alphas and update score and importance for 
    unsigned int nAlphas = static_cast<unsigned int>(alphaMinVector.size());

    for(unsigned int i = 0u; i < nAlphas; i++) {

      double currAlpha{static_cast<double>(i)/(nAlphas-1)};

      for(auto tvinfo : tempVertexInfoMap) {

        double ltc,mtc,ltg,mtg;
        Vertex lpred,mpred,lsucc,msucc;
        std::tie(ltc,lpred,mtc,mpred,ltg,lsucc,mtg,msucc) = tvinfo.second;

        double alphaScore{currAlpha*(ltc+ltg) + (1.0 - currAlpha)*(mtc+mtg)};

        // Subtract importance if it gets worse
        
        if( alphaScore < alphaMinVector[i]) {
          tempAlphaVector[i] = alphaScore;
          tempAlphaIdxVector[i] = tvinfo.first;
          tempVertexImportanceMap[tvinfo.first] += alphaMinVector[i] - alphaScore;
        }
        else if (alphaScore < alphaSecondMinVector[i]) {
          tempSecondAlphaVector[i] = alphaScore;
        }
      }
    }

    double newRoadmapScore{getScoreFromAlphaVector(tempAlphaVector,tempAlphaIdxVector)};
    double improvement{0.0};

    // Check if roadmap score improves
    if(mCurrRoadmapScore - newRoadmapScore > 0.00001) {

      mSuccPerturbations++;
      improvement = mCurrRoadmapScore - newRoadmapScore;
      mCurrRoadmapScore = newRoadmapScore;

      // Change state of vertex and other values
      mSpace->copyState(mCurrentRoadmap[u].v_state->state, _perturbedState->state);

      L_to_come[u] = min_L_tc;
      M_to_come[u] = min_M_tc;
      L_to_go[u] = min_L_tg;
      M_to_go[u] = min_M_tg;
      L_preds[u] = L_pred;
      M_preds[u] = M_pred;
      L_succs[u] = L_succ;
      M_succs[u] = M_succ;


      // Update edges
      for (auto nbr_wts : newNbrVertexWeights)
      {
        Vertex nbr{nbr_wts.first};
        std::pair<double,double> weights{nbr_wts.second};

        // If Edge exists, either update or leave it if weight inf
        std::pair<Edge,bool> potential_edge = boost::edge(u,nbr,mCurrentRoadmap);

        if(potential_edge.second) {
          mCurrentRoadmap[potential_edge.first].distance = weights.first;
          mCurrentRoadmap[potential_edge.first].collMeasure = weights.second;

          if(weights.first == std::numeric_limits<double>::infinity()) {
            mInfiniteCostEdges.insert(potential_edge.first);
          }
        }
        else{
          std::pair<Edge,bool> new_edge = boost::add_edge(u,nbr,mCurrentRoadmap);
          mCurrentRoadmap[new_edge.first].distance = weights.first;
          mCurrentRoadmap[new_edge.first].collMeasure = weights.second;
          mCurrentRoadmap[new_edge.first].blockedStatus = UNKNOWN;
          mCurrentRoadmap[new_edge.first].hasPoints = false;
        }
      }


      // Update other vertices
      for(auto tvinfo : tempVertexInfoMap)
      {
        double ltc,mtc,ltg,mtg;
        Vertex lpred,mpred,lsucc,msucc;
        std::tie(ltc,lpred,mtc,mpred,ltg,lsucc,mtg,msucc) = tvinfo.second;

        L_to_come[tvinfo.first] = ltc;
        M_to_come[tvinfo.first] = mtc;
        L_to_go[tvinfo.first] = ltg;
        M_to_go[tvinfo.first] = mtg;
        L_preds[tvinfo.first] = lpred;
        M_preds[tvinfo.first] = mpred;
        L_succs[tvinfo.first] = lsucc;
        M_succs[tvinfo.first] = msucc;

        vertexImportance[tvinfo.first] = tempVertexImportanceMap[tvinfo.first];
      }


      // Update alpha min vector and second min vector
      alphaMinVector = tempAlphaVector;
      alphaMinIdxVector = tempAlphaIdxVector;
      alphaSecondMinVector = tempSecondAlphaVector;
    }
    
    // Either way, remove new vertex
    boost::remove_vertex(tempNewVertex, mCurrentRoadmap);

    // Return numerical value of improvement
    return improvement;

  }


  // NOTE - Just assuming start and goal are 0/1 or 1/0 respectively
  unsigned int softMaxVertexSelection()
  {
    // Uses the current alpha score vector to select a vertex
    std::vector<double> softmaxScores(mNumCurrVertices - 2);
    double impScoreSum{0.0};

    for(unsigned int i=2u ; i < mNumCurrVertices; i++)
    {
      impScoreSum += vertexImportance[i];
    }

    double softmaxScoreSum{0.0};
    for(unsigned int i=2u ; i < mNumCurrVertices; i++)
    {
      softmaxScores[i-2] = std::exp(-vertexImportance[i]/impScoreSum);
      softmaxScoreSum += softmaxScores[i-2];
    }

    std::uniform_real_distribution<> dis(0,1);
    double threshold{dis(mGen)};
    double thresholdCheck{0.0};

    for(unsigned int idx=2u ; idx < mNumCurrVertices; idx++)
    {
      thresholdCheck += softmaxScores[idx-2]/softmaxScoreSum;

      if(thresholdCheck > threshold) {
        return idx;
      }
    }

    return (mNumCurrVertices - 1);
  }


  unsigned int randomVertexSelection()
  {
    int numVerts = static_cast<int>(mNumCurrVertices);
    std::uniform_int_distribution<> dis(2,mNumCurrVertices-1);
    return static_cast<unsigned int>(dis(mGen));
  }


// TODO : Remove all edges with either weight infinite in main loop

  void updateRoadmap()
  {
    // First add initial samples
    addInitialBatch();

    // If simple flag, just return
    if(mSimpleFlag) {
      return;
    }

    // Compute initial score
    computeInitialScore();

    for(unsigned int i=0; i < mNumPerturbations; i++)
    {
      // Choose a sample based on its importance
      Vertex chosenVert{randomVertexSelection()};
      StateConPtr perturbedState(std::make_shared<StateCon>(mSpace));
      if(perturbVertexApproxGrad(chosenVert,perturbedState)){
        double improvement = implementPerturbation(chosenVert, perturbedState);
        if(improvement > 0.00001){
          std::cout<<"Trial "<<i<<" : Vertex "<<chosenVert<<" improved score by "<<improvement<<std::endl;
        }
      }
    }

    std::cout<<"Final score - "<<mCurrRoadmapScore<<std::endl;
    std::cout<<mSuccPerturbations<<" successful perturbations"<<std::endl;

    // TODO : See if you need this
    for(Edge e : mInfiniteCostEdges)
    {
      if(mCurrentRoadmap[e].distance == std::numeric_limits<double>::infinity()) {
        boost::remove_edge(e,mCurrentRoadmap);
      }
    }
  
  }


private:

  // Testing - step size of perturbation
  double mPerturbSize;
  unsigned int mNumPerturbations;

  // DEBUGGING
  unsigned int mSuccPerturbations;

  std::mt19937 mGen;

  bool mSimpleFlag;

  Graph& mFullRoadmap;
  ompl::NearestNeighbors<Vertex>& mCurrVertexNN;

  VertexIter mCurrVertex;
  VertexIter mLastVertex;
  Vertex mStartVertex;
  Vertex mGoalVertex;
  unsigned int mNumCurrVertices;

  double mProbThreshold;

  std::set<Edge> mInfiniteCostEdges; 

  // For tracking additional info for vertices
  std::map<Vertex,double> M_to_come;
  std::map<Vertex,Vertex> M_preds;
  std::map<Vertex,double> L_to_come;
  std::map<Vertex,Vertex> L_preds;
  std::map<Vertex,double> M_to_go;
  std::map<Vertex,Vertex> M_succs;
  std::map<Vertex,double> L_to_go;
  std::map<Vertex,Vertex> L_succs;
  std::vector<double> vertexImportance;
  std::vector<double> alphaMinVector;
  std::vector<unsigned int> alphaMinIdxVector;
  std::vector<double> alphaSecondMinVector;

  double mCurrRoadmapScore;


};


}
}

#endif
