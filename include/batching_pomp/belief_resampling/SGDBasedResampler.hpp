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
#include <algorithm>

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

static const double STEPSIZE_EPSILON = 0.001;
static const double GRAD_EPSILON = 0.001;

class EdgeSearchWeightMap
{
public:
  typedef boost::readable_property_map_tag category;
  typedef Edge key_type;
  typedef double value_type;
  typedef double reference;
  Graph& mCurrRoadmap;
  double mAlpha;
  EdgeSearchWeightMap(Graph& _roadmap, double _alpha)
  : mCurrRoadmap(_roadmap), mAlpha{_alpha} {}
};
const double get(const EdgeSearchWeightMap& _ewMap, const Edge& e)
{
  double w_m{_ewMap.mCurrRoadmap[e].collMeasure};
  double w_l{_ewMap.mCurrRoadmap[e].distance};

  double exp_cost{_ewMap.mAlpha*w_l + (1-_ewMap.mAlpha)*w_m};
  return exp_cost;
}

class throw_visitor_general
{

public:
  Vertex mVThrow;
  throw_visitor_general(Vertex _vThrow)
  : mVThrow{_vThrow} 
  {}
  inline void initialize_vertex(Vertex u, const Graph& g) {}
  inline void discover_vertex(Vertex u, const Graph& g) {}
  inline void examine_vertex(Vertex u, const Graph& g)
  {
    if(u == mVThrow) {
      throw throw_visitor_exception();
    }
  }
  inline void examine_edge(Edge e, const Graph& g){}
  inline void edge_relaxed(Edge e, const Graph & g) {}
  inline void edge_not_relaxed(Edge e, const Graph & g) {}
  inline void black_target(Edge e, const Graph & g) {}
  inline void finish_vertex(Vertex u, const Graph & g) {}
};



class SGDBasedResampler : public virtual BeliefInformedResampler
{

public:

  SGDBasedResampler(BatchingPOMP& _planner,
                    unsigned int _numPerturbations,
                    double _perturbSize,
                    double _probThreshold = 1.0,
                    double _randomOffsetSize=0.0)
  : BeliefInformedResampler(_planner)
  , mFullRoadmap(*(_planner.full_g))
  , mCurrVertexNN(*(_planner.mVertexNN))
  , mStartVertex{_planner.getStartVertex()}
  , mGoalVertex{_planner.getGoalVertex()}
  , mNumPerturbations{_numPerturbations}
  , mPerturbSize{_perturbSize}
  , mProbThreshold{_probThreshold}
  , mSuccPerturbations{0u}
  , mRandomOffsetSize{_randomOffsetSize}
  {
    boost::tie(mCurrVertex,mLastVertex) = vertices(mFullRoadmap);
    std::random_device rd;
    mGen.seed(rd());
  }

  /// Add the next n samples from list subject to thresholding conditions
  void addInitialBatch()
  {
    // Iterate and add to current roadmap
    unsigned int numVerticesAdded{0u};
    std::vector<Vertex> vertex_vector(mBatchParams.first);

    while(numVerticesAdded < mBatchParams.first)
    {
      StateConPtr perturbedState(std::make_shared<StateCon>(mSpace));
      if(mRandomOffsetSize > 0.0){
        ompl::base::RealVectorStateSampler rvSampler(mSpace.get());
        rvSampler.sampleUniformNear(perturbedState->state, mFullRoadmap[*mCurrVertex].v_state->state, mRandomOffsetSize);
      }
      else{
        mSpace->copyState(perturbedState->state,mFullRoadmap[*mCurrVertex].v_state->state);
      }

      cspacebelief::BeliefPoint query(perturbedState->state, 
                        mSpace->getDimension(), -1.0);

      if(mPlanner.getSpaceInformation()->isValid(perturbedState->state) &&
         mBeliefModel->estimate(query) < mProbThreshold) {
        Vertex newVertex{boost::add_vertex(mCurrentRoadmap)};

        mCurrentRoadmap[newVertex].v_state = perturbedState;
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

    mVertexImportance.resize(mNumCurrVertices,0.0);

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

  double getCurrRoadmapScore() const
  {
    return mCurrRoadmapScore;
  }

  double getParetoConvexHullScore(std::vector<double>& _tempVertexImportance,
                                  std::map<Vertex,std::set<std::pair<double,double>>>& _tempVertexToAlphaRangeMap)
  {
    // ASSUMES tempVertexImportance initialized with zeros

    typedef std::tuple<double,std::deque<Vertex>,std::pair<double, double> > alphaPathWeightTuple;

    std::set<alphaPathWeightTuple> uniqueCvxHullPaths;

    std::deque< std::pair<alphaPathWeightTuple, alphaPathWeightTuple> > frontierQ;

    std::vector< alphaPathWeightTuple > zeroOnePaths(2);


    std::function<double(const Vertex& u, const Vertex& v)> vertDistFn = 
      std::bind(&SGDBasedResampler::vertexDistFun,this,std::placeholders::_1,std::placeholders::_2);

    //First try for alpha 0 and 1
    for(unsigned int i = 0u; i <= 1u; i++)
    {
      double currAlpha = static_cast<double>(i);
      std::unique_ptr<boost::astar_heuristic<Graph, double>> heuristicFn;
      if(i == 0u){
        heuristicFn.reset(new zero_heuristic<Graph,double>());
      }
      else{
        heuristicFn.reset(new exp_distance_heuristic<Graph,double>(1.0,mGoalVertex,vertDistFn));
      }

      std::map<Vertex,Vertex> startPreds;
      std::map<Vertex,double> startDist;
      std::map<Vertex,double> startFValue;
      std::map<Vertex,boost::default_color_type> colorMap;

      try
      {
        boost::astar_search(
          mCurrentRoadmap,
          mStartVertex,
          *heuristicFn,
          throw_visitor_general(mGoalVertex),
          boost::make_assoc_property_map(startPreds),
          boost::make_assoc_property_map(startFValue),
          boost::make_assoc_property_map(startDist),
          EdgeSearchWeightMap(mCurrentRoadmap,currAlpha),
          boost::get(boost::vertex_index, mCurrentRoadmap),
          boost::make_assoc_property_map(colorMap),
          std::less<double>(),
          boost::closed_plus<double>(std::numeric_limits<double>::max()),
          std::numeric_limits<double>::max(),
          double()
        );
      }
      catch (const throw_visitor_exception & ex)
      {
      }

      if(startDist[mGoalVertex] == std::numeric_limits<double>::infinity()) {
        throw ompl::Exception("No candidate path found in roadmap");
      }

      std::deque<Vertex> thisAlphaPath;

      std::pair<double, double> pathWeights = 
          getPathWeightsAfterSearch(startPreds, thisAlphaPath);
      
      zeroOnePaths[i] = std::make_tuple(currAlpha, thisAlphaPath, pathWeights);
    }

    // Check if zero-one alpha equal
    if(std::get<1>(zeroOnePaths[0]) == std::get<1>(zeroOnePaths[1]))
    {
      uniqueCvxHullPaths.insert(zeroOnePaths[0]);
    }
    else
    {
      // Now do the binary splitting stuff
      uniqueCvxHullPaths.insert(zeroOnePaths[0]);
      uniqueCvxHullPaths.insert(zeroOnePaths[1]);
      frontierQ.push_back( std::make_pair(zeroOnePaths[0], zeroOnePaths[1]) );

      while(!frontierQ.empty())
      {
        std::pair<alphaPathWeightTuple, alphaPathWeightTuple> 
            frontierQtop = frontierQ.front();

        frontierQ.pop_front();

        double alphaLow, alphaHigh;
        std::deque<Vertex> pathLow, pathHigh;
        std::pair<double,double> weightsLow, weightsHigh;

        std::tie(alphaLow,pathLow,weightsLow) = frontierQtop.first;
        std::tie(alphaHigh,pathHigh,weightsHigh) = frontierQtop.second;

        double newCurrAlpha = 
          (weightsHigh.second - weightsLow.second)/(weightsLow.first - weightsHigh.first + weightsHigh.second - weightsLow.second);


        // Now search with newAlpha
        std::map<Vertex,Vertex> startPreds;
        std::map<Vertex,double> startDist;
        std::map<Vertex,double> startFValue;
        std::map<Vertex,boost::default_color_type> colorMap;

        std::unique_ptr<boost::astar_heuristic<Graph, double>> heuristicFn;
        heuristicFn.reset(new exp_distance_heuristic<Graph,double>(newCurrAlpha,mGoalVertex,vertDistFn));

        try
        {
          boost::astar_search(
            mCurrentRoadmap,
            mStartVertex,
            *heuristicFn,
            throw_visitor_general(mGoalVertex),
            boost::make_assoc_property_map(startPreds),
            boost::make_assoc_property_map(startFValue),
            boost::make_assoc_property_map(startDist),
            EdgeSearchWeightMap(mCurrentRoadmap,newCurrAlpha),
            boost::get(boost::vertex_index, mCurrentRoadmap),
            boost::make_assoc_property_map(colorMap),
            std::less<double>(),
            boost::closed_plus<double>(std::numeric_limits<double>::max()),
            std::numeric_limits<double>::max(),
            double()
          );
        }
        catch (const throw_visitor_exception & ex)
        {
        }

        if(startDist[mGoalVertex] == std::numeric_limits<double>::max()) {
          throw ompl::Exception("No candidate path found in roadmap");
        }

        std::deque<Vertex> newCurrAlphaPathVertices;
        std::pair<double, double> newCurrAlphaPathWeights = 
          getPathWeightsAfterSearch(startPreds, newCurrAlphaPathVertices);

        // Now check if different from previous
        if(newCurrAlphaPathVertices == pathLow || newCurrAlphaPathVertices == pathHigh)
        {
          // Same, break
          break;
        }
        else{
          // Now enter new path in set and update frontier
          alphaPathWeightTuple newFrontierElem = 
            std::make_tuple(newCurrAlpha, newCurrAlphaPathVertices, newCurrAlphaPathWeights);
          uniqueCvxHullPaths.insert(newFrontierElem);

          //Now create two frontier elements with low and new and new and hi
          frontierQ.push_back( std::make_pair(frontierQtop.first, newFrontierElem));
          frontierQ.push_back( std::make_pair(newFrontierElem, frontierQtop.second));
        }
      }
    }

    double currAlpha{0.0};
    double nextAlpha{0.0};
    auto it_next = uniqueCvxHullPaths.begin();
    it_next++;
    auto it = uniqueCvxHullPaths.begin();

    double roadmapScore{0.0};

    //Now determine what to do with uniqueCvxHullPaths
    for(; it_next != uniqueCvxHullPaths.end(); it++, it_next++)
    {
      std::pair<double,double> currWeights = std::get<2>(*it);
      std::pair<double,double> nextWeights = std::get<2>(*it_next);

      // Now compute the next alpha with atan2
      nextAlpha = (nextWeights.second - currWeights.second)/(currWeights.first - nextWeights.first + nextWeights.second - currWeights.second);


      // So this path is important from currAlpha to nextAlpha
      std::deque<Vertex> currAlphaPathVertices = std::get<1>(*it);

      double pathContrib = 
        (currWeights.first/2.0)*(nextAlpha*nextAlpha - currAlpha*currAlpha) + 
        (currWeights.second)*( (nextAlpha - 0.5*nextAlpha*nextAlpha) - (currAlpha - 0.5*currAlpha*currAlpha) );

      if(pathContrib < std::numeric_limits<double>::epsilon()){
        continue;
      }

      roadmapScore += pathContrib;

      //std::cout<<"Pathcontrib "<<pathContrib<<" from "<<currAlpha<<" to "<<nextAlpha<<std::endl;

      for(Vertex pv : currAlphaPathVertices)
      {
        //std::cout<<pv<<" ";
        _tempVertexImportance[pv] += 1.0/pathContrib;
        auto search = _tempVertexToAlphaRangeMap.find(pv);
        if(search != _tempVertexToAlphaRangeMap.end()){
          search->second.insert( std::make_pair(currAlpha, nextAlpha));
        }
        else{
          _tempVertexToAlphaRangeMap.insert(std::make_pair(pv, std::set<std::pair<double,double>>()));
          _tempVertexToAlphaRangeMap[pv].insert( std::make_pair(currAlpha, nextAlpha));
        }
      }

      currAlpha = nextAlpha;
      //std::cout<<"Weights : "<<currWeights.first<<","<<currWeights.second<<std::endl;
    }

    // Now do for the alpha = 1 case
    nextAlpha = 1.0;
    std::pair<double,double> lastWeights = std::get<2>(*it);
    std::deque<Vertex> lastAlphaPathVertices = std::get<1>(*it);
    double pathContrib = 
        (lastWeights.first/2.0)*(nextAlpha*nextAlpha - currAlpha*currAlpha) + 
        (lastWeights.second)*( (nextAlpha - 0.5*nextAlpha*nextAlpha) - (currAlpha - 0.5*currAlpha*currAlpha) );
    roadmapScore += pathContrib;
    //std::cout<<"Pathcontrib "<<pathContrib<<" from "<<currAlpha<<" to "<<nextAlpha<<" ; path - ";
    for(Vertex pv : lastAlphaPathVertices)
    {
      //std::cout<<pv<<" ";
      _tempVertexImportance[pv] += 1.0/pathContrib;
      auto search = _tempVertexToAlphaRangeMap.find(pv);
      if(search != _tempVertexToAlphaRangeMap.end()){
        search->second.insert( std::make_pair(currAlpha, nextAlpha));
      }
      else{
        _tempVertexToAlphaRangeMap.insert(std::make_pair(pv, std::set<std::pair<double,double>>()));
        _tempVertexToAlphaRangeMap[pv].insert( std::make_pair(currAlpha, nextAlpha));
      }
    }
    //std::cout<<"Weights : "<<lastWeights.first<<","<<lastWeights.second<<std::endl;
    return roadmapScore;

  }

  double vertexDistFun(const Vertex& u, const Vertex& v) const
  {
    return mSpace->distance(mCurrentRoadmap[u].v_state->state, mCurrentRoadmap[v].v_state->state);
  }

  void computeInitialScore()
  {
    mCurrRoadmapScore = getParetoConvexHullScore(mVertexImportance, mVertexToAlphaRangeMap);
  }

  std::pair<double,double> getPathWeightsAfterSearch(const std::map<Vertex, Vertex>& _startPreds,
                                                     std::deque<Vertex>& _pathVertices)
  {
    // We are assuming that a path exists from start to goal
    Vertex vWalk{mGoalVertex};
    _pathVertices.push_front(mGoalVertex);
    double pathLength{0.0};
    double pathMeasure{0.0};

    while(vWalk!=mStartVertex)
    {
      Vertex vPred{_startPreds.at(vWalk)};
      _pathVertices.push_front(vPred);
      std::pair<Edge, bool> edgePair{boost::edge(vPred, vWalk, mCurrentRoadmap)};
      if(!edgePair.second) {
        throw ompl::Exception("Error! Edge present during search no longer exists in graph!");
      }

      pathLength += mCurrentRoadmap[edgePair.first].distance;
      pathMeasure += mCurrentRoadmap[edgePair.first].collMeasure;
      vWalk = vPred;
    }

    return std::make_pair(pathLength, pathMeasure);
  }

  double getScoreFromAlphaVector(const std::vector<std::pair<double,double>>& _alphaWeightsVector) const
  {
    double score{_alphaWeightsVector[0].second};
    double prevMeasure{_alphaWeightsVector[0].second};
    double prevLength{_alphaWeightsVector[0].first};
    unsigned int nAlphas = static_cast<unsigned int>(_alphaWeightsVector.size());
    unsigned int numUniqueVals{1u};

    for(unsigned int i = 1u; i < nAlphas; i++) {

      double currAlpha{static_cast<double>(i)/(nAlphas-1)};
      score += currAlpha*_alphaWeightsVector[i].first + (1.0 - currAlpha)*_alphaWeightsVector[i].second;
      if(std::fabs(_alphaWeightsVector[i].first - prevLength) > std::numeric_limits<double>::epsilon() 
        || std::fabs(_alphaWeightsVector[i].second - prevMeasure) > std::numeric_limits<double>::epsilon())
      {
        numUniqueVals ++;
        prevLength = _alphaWeightsVector[i].first;
        prevMeasure = _alphaWeightsVector[i].second;
      }
    }

    // TODO : Play around with this
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

  bool perturbVertexApproxGradImportance(const Vertex& u, StateConPtr& perturbedState) const
  {
    double alphaWeightSum{0.0};
    
    std::map<Edge,std::set<std::pair<double, double>>> importantEdgesWithAlphaRanges;

    // First find all neighbours who have non-zero importance and store those edges
    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei) {
      Vertex nbr{target(*ei,mCurrentRoadmap)};
      if(mVertexImportance[target(*ei,mCurrentRoadmap)] > 0){
        
        // Find intersection between alpha sets of u and nbr 
        // And add those to edgeAlphas
        std::set<std::pair<double, double>> alpha_range_intersection;
        std::set_intersection(mVertexToAlphaRangeMap.at(u).begin(), mVertexToAlphaRangeMap.at(u).end(),
                              mVertexToAlphaRangeMap.at(nbr).begin(), mVertexToAlphaRangeMap.at(nbr).end(),
                              std::inserter(alpha_range_intersection,alpha_range_intersection.begin()));
        if(alpha_range_intersection.size() > 0){
          importantEdgesWithAlphaRanges.insert(std::make_pair(*ei,alpha_range_intersection));
        }
      }
    }

    // Iterate over alpha indices which vertex minimized
    //std::cout<<"Alpha Ranges to care about - ";
    for(auto edgeAlphaRanges : importantEdgesWithAlphaRanges)
    {
      Edge e{edgeAlphaRanges.first};
      std::set<std::pair<double, double>> impAlphaRanges(edgeAlphaRanges.second);
      //std::cout<<e<<" : ("<<mCurrentRoadmap[e].distance<<","<<mCurrentRoadmap[e].collMeasure<<") - ";

      for(auto alphaRange : impAlphaRanges) {
        //std::cout<<alphaRange.first<<" to "<<alphaRange.second<<"  ";
        alphaWeightSum +=
          (mCurrentRoadmap[e].distance/2.0)*(alphaRange.second*alphaRange.second - alphaRange.first*alphaRange.first) +
          (mCurrentRoadmap[e].collMeasure)*( (alphaRange.second - 0.5*alphaRange.second*alphaRange.second) 
                                           - (alphaRange.first - 0.5*alphaRange.first*alphaRange.first) );
      }
      //std::cout<<std::endl;
    }

    // Get estimate of vertex
    unsigned int nDims{mSpace->getDimension()};
    cspacebelief::BeliefPoint query(mCurrentRoadmap[u].v_state->state, 
                                    nDims, -1.0);
    double vertLogProb{-std::log(1.0 - mBeliefModel->estimate(query))};
    Eigen::VectorXd grads(nDims);

    //std::cout<<"Alphawtsum vs new wts : "<<alphaWeightSum<<std::endl;

    bool localMin{true};

    for(unsigned int d=0u; d < nDims; d++)
    {
      StateConPtr oneDimPerturbed(std::make_shared<StateCon>(mSpace));
      mSpace->copyState(oneDimPerturbed->state, mCurrentRoadmap[u].v_state->state);
      // Perturb along dimension
      double* values = oneDimPerturbed->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      values[d] = values[d] + GRAD_EPSILON;

      StateConPtr oneDimPerturbed2(std::make_shared<StateCon>(mSpace));
      mSpace->copyState(oneDimPerturbed2->state, mCurrentRoadmap[u].v_state->state);
      // Perturb along dimension
      double* values2 = oneDimPerturbed2->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      values2[d] = values2[d] - GRAD_EPSILON;

      double newAlphaWtSum{0.0};
      double newAlphaWtSum2{0.0};

      for(auto edgeAlphaRanges : importantEdgesWithAlphaRanges)
      {
        Edge e{edgeAlphaRanges.first};
        std::set<std::pair<double, double>> impAlphaRanges(edgeAlphaRanges.second);

        Vertex v{boost::target(e,mCurrentRoadmap)};

        double new_length{mSpace->distance(oneDimPerturbed->state, mCurrentRoadmap[v].v_state->state)};
        cspacebelief::BeliefPoint per_query(oneDimPerturbed->state, 
                      mSpace->getDimension(), -1.0);
        double newLogProb{-std::log(1.0 - mBeliefModel->estimate(per_query))};
        double new_meas = mCurrentRoadmap[e].collMeasure + newLogProb - vertLogProb;

        double new_length2{mSpace->distance(oneDimPerturbed2->state, mCurrentRoadmap[v].v_state->state)};
        cspacebelief::BeliefPoint per_query2(oneDimPerturbed2->state, 
                      mSpace->getDimension(), -1.0);
        double newLogProb2{-std::log(1.0 - mBeliefModel->estimate(per_query2))};
        double new_meas2 = mCurrentRoadmap[e].collMeasure + newLogProb2 - vertLogProb;

        // std::cout<<e<<" : ("<<new_length<<","<<new_meas<<") ";
        // std::cout<<e<<" : ("<<new_length2<<","<<new_meas2<<") ";

        for(auto alphaRange : impAlphaRanges) {
          //std::cout<<alphaRange.first<<" to "<<alphaRange.second;
          newAlphaWtSum += 
            (new_length/2.0)*(alphaRange.second*alphaRange.second - alphaRange.first*alphaRange.first) +
            (new_meas)*( (alphaRange.second - 0.5*alphaRange.second*alphaRange.second) 
                                           - (alphaRange.first - 0.5*alphaRange.first*alphaRange.first) );
          newAlphaWtSum2 += 
            (new_length2/2.0)*(alphaRange.second*alphaRange.second - alphaRange.first*alphaRange.first) +
            (new_meas2)*( (alphaRange.second - 0.5*alphaRange.second*alphaRange.second) 
                                           - (alphaRange.first - 0.5*alphaRange.first*alphaRange.first) );
        }
      //std::cout<<std::endl;  
      }

      //std::cout<<newAlphaWtSum<<" ; "<<newAlphaWtSum2<<std::endl;
      grads[d] = (newAlphaWtSum - newAlphaWtSum2) / (2*GRAD_EPSILON);

      if(newAlphaWtSum < alphaWeightSum || newAlphaWtSum2 < alphaWeightSum){
        localMin = false;
      }
    }

    if(localMin){
      return false;
    }

    //std::cout<<std::endl<<grads.norm()<<std::endl;
    grads.normalize();

    if(std::isnan(grads[0]) || std::isnan(grads[1])){
      return false;
    }

    //std::cout<<grads<<std::endl;

    //throw ompl::Exception("IMP STOP!");

    double stepSize{mPerturbSize};
    StateConPtr perturbedStateCopy(std::make_shared<StateCon>(mSpace));
    while(stepSize > STEPSIZE_EPSILON)
    {
      mSpace->copyState(perturbedStateCopy->state, mCurrentRoadmap[u].v_state->state);
      double* values = perturbedStateCopy->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      for(unsigned int d=0u; d < nDims; d++)
      {
        values[d] = values[d] - stepSize*grads[d];
      }

      if(mSpace->satisfiesBounds(perturbedStateCopy->state) &&
            mPlanner.getSpaceInformation()->isValid(perturbedStateCopy->state)) {
        break;
      }

      stepSize /= 2.0;
    }
    if(stepSize > STEPSIZE_EPSILON) {
      mSpace->copyState(perturbedState->state,perturbedStateCopy->state);
      return true;
    }

    return false;

  }

  bool perturbVertexApproxGradNoImportance(const Vertex& u, StateConPtr& perturbedState) const
  {
    double edgeWeightSum{0.0};

    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei)
    {
      double length{mCurrentRoadmap[*ei].distance};
      double meas{mCurrentRoadmap[*ei].collMeasure};
      
      edgeWeightSum += length + meas;
    }

    // Get estimate of vertex
    unsigned int nDims{mSpace->getDimension()};
    cspacebelief::BeliefPoint query(mCurrentRoadmap[u].v_state->state, 
                                    nDims, -1.0);
    double vertLogProb{-std::log(1.0 - mBeliefModel->estimate(query))};
    Eigen::VectorXd grads(nDims);

    bool localMin{true};

    for(unsigned int d=0u; d < nDims; d++)
    {
      StateConPtr oneDimPerturbed(std::make_shared<StateCon>(mSpace));
      mSpace->copyState(oneDimPerturbed->state, mCurrentRoadmap[u].v_state->state);

      // Perturb along dimension
      double* values = oneDimPerturbed->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      values[d] = values[d] + GRAD_EPSILON;

      StateConPtr oneDimPerturbed2(std::make_shared<StateCon>(mSpace));
      mSpace->copyState(oneDimPerturbed2->state, mCurrentRoadmap[u].v_state->state);
      // Perturb along dimension
      double* values2 = oneDimPerturbed2->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      values2[d] = values2[d] - GRAD_EPSILON;

      double newEdgeWtSum{0.0};
      double newEdgeWtSum2{0.0};

      for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei)
      {
        Vertex v{boost::target(*ei,mCurrentRoadmap)};

        double new_length{mSpace->distance(oneDimPerturbed->state, mCurrentRoadmap[v].v_state->state)};
        cspacebelief::BeliefPoint per_query(oneDimPerturbed->state, 
                        mSpace->getDimension(), -1.0);
        double newLogProb{-std::log(1.0 - mBeliefModel->estimate(per_query))};
        double new_meas = mCurrentRoadmap[*ei].collMeasure + newLogProb - vertLogProb;
        newEdgeWtSum += new_length + new_meas;

        double new_length2{mSpace->distance(oneDimPerturbed2->state, mCurrentRoadmap[v].v_state->state)};
        cspacebelief::BeliefPoint per_query2(oneDimPerturbed2->state, 
                      mSpace->getDimension(), -1.0);
        double newLogProb2{-std::log(1.0 - mBeliefModel->estimate(per_query2))};
        double new_meas2 = mCurrentRoadmap[*ei].collMeasure + newLogProb2 - vertLogProb;

        newEdgeWtSum2 += new_length2 + new_meas2;

      }

      grads[d] = (newEdgeWtSum - newEdgeWtSum2) / (2.0*GRAD_EPSILON);

      if(newEdgeWtSum < edgeWeightSum || newEdgeWtSum2 < edgeWeightSum){
        localMin = false;
      }

    }

    if(localMin){
      return false;
    }

    grads.normalize();
    if(std::isnan(grads[0]) || std::isnan(grads[1])){
      return false;
    }

    //std::cout<<"For vertex of no importance - grads : "<<grads<<std::endl;

    //throw ompl::Exception("NONIMP STOP!");

    double stepSize{mPerturbSize};
    StateConPtr perturbedStateCopy(std::make_shared<StateCon>(mSpace));
    while(stepSize > STEPSIZE_EPSILON)
    {
      mSpace->copyState(perturbedStateCopy->state, mCurrentRoadmap[u].v_state->state);
      double* values = perturbedStateCopy->state->as<
        ompl::base::RealVectorStateSpace::StateType>()->values;
      for(unsigned int d=0u; d < nDims; d++)
      {
        values[d] = values[d] - stepSize*grads[d];
      }

      if( mSpace->satisfiesBounds(perturbedStateCopy->state) &&
            mPlanner.getSpaceInformation()->isValid(perturbedStateCopy->state)) {
        break;
      }

      stepSize /= 2.0;
    }

    if(stepSize > STEPSIZE_EPSILON) {
      mSpace->copyState(perturbedState->state,perturbedStateCopy->state);
      return true;
    }

    return false;
  }

  double implementPerturbation(Vertex& u, const StateConPtr& _perturbedState)
  {
    // TODO - Put back threshold check? No, do it in approx grad
    // ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
    //   old_state(mSpace, mCurrentRoadmap[u].v_state->state);
    // std::cout<<"Old state - "<<old_state<<std::endl;
    // ompl::base::ScopedState<ompl::base::RealVectorStateSpace>
    //   new_state(mSpace, _perturbedState->state);
    // std::cout<<"New state - "<<new_state<<std::endl;


    // First store old state 
    StateConPtr oldState = std::make_shared<StateCon>(mSpace);
    mSpace->copyState(oldState->state, mCurrentRoadmap[u].v_state->state);

    // Now change vertex and recompute lengths and measures
    std::map<Vertex,std::pair<double,double>> oldEdgeMap;
    mSpace->copyState(mCurrentRoadmap[u].v_state->state, _perturbedState->state);

    std::set<Edge> oldEdgesToRemove;

    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,mCurrentRoadmap); ei!=ei_end; ++ei)
    {
      Vertex nbr{boost::target(*ei,mCurrentRoadmap)};
      // Compute new weights and store old weights
      oldEdgeMap.insert(std::make_pair(nbr,
        std::make_pair(mCurrentRoadmap[*ei].distance, mCurrentRoadmap[*ei].collMeasure)));

      double edgeDist{vertexDistFun(u,nbr)};

      if(edgeDist > mBatchParams.second)
      {
        mCurrentRoadmap[*ei].distance = std::numeric_limits<double>::infinity();
        mCurrentRoadmap[*ei].collMeasure = std::numeric_limits<double>::infinity();
        oldEdgesToRemove.insert(*ei);
      }
      else{
        mCurrentRoadmap[*ei].distance = edgeDist;
        mCurrentRoadmap[*ei].collMeasure = mPlanner.computeEdgeCollisionMeasureNoStates(u,nbr);
      }
    }

    // Now add potential new nbrs with their weights
    // TODO - Check that this actually returns new neighbours!
    std::set<Edge> newEdgesToRemove;
    std::vector<Vertex> tempVertexNbrs;
    mCurrVertexNN.nearestR(u,mBatchParams.second,tempVertexNbrs);

    for(Vertex tvnbr : tempVertexNbrs) {

      if(tvnbr != u && !boost::edge(u,tvnbr,mCurrentRoadmap).second) {
        //TODO : Check this happens at least once for big perturbations
        std::pair<Edge,bool> new_edge = boost::add_edge(u,tvnbr,mCurrentRoadmap);
        mCurrentRoadmap[new_edge.first].distance = mPlanner.vertexDistFun(u,tvnbr);
        mCurrentRoadmap[new_edge.first].collMeasure = 
          mPlanner.computeEdgeCollisionMeasureNoStates(u,tvnbr);
        mCurrentRoadmap[new_edge.first].blockedStatus = UNKNOWN;
        mCurrentRoadmap[new_edge.first].hasPoints = false;
        
        newEdgesToRemove.insert(new_edge.first);
      }
    }


    // Bookkeeping stuff
    std::vector<double> tempVertexImportance(mNumCurrVertices,0.0);
    std::map<Vertex,std::set<std::pair<double,double>>> tempVertexToAlphaRangeMap;

    double newRoadmapScore{getParetoConvexHullScore(tempVertexImportance,tempVertexToAlphaRangeMap)};

    double improvement{0.0};

    // if(mVertexImportance[u] > 0.0){
    //   std::cout<<"For vertex "<<u<<" of importance "<<mVertexImportance[u]<<", the improvement is "
    //   <<(mCurrRoadmapScore-newRoadmapScore)<<std::endl;
    //   //throw ompl::Exception("STOP!");
    // }
      
    // Now compare with previous and accept if improved
    if(mCurrRoadmapScore >= newRoadmapScore)
    {
      improvement = mCurrRoadmapScore - newRoadmapScore;
      mSuccPerturbations++;
      // And copy weight vectors and maps
      mCurrRoadmapScore = newRoadmapScore;
      mVertexImportance = tempVertexImportance;
      mVertexToAlphaRangeMap = tempVertexToAlphaRangeMap;
      
      // Either way implement change
      // Which is just to remove oldEdgesToRemove
      for(Edge e : oldEdgesToRemove){
        boost::remove_edge(e,mCurrentRoadmap);
      }
    }
    else
    {
      // We ended up reducing the score - don't accept changes
      // Restore state
      mSpace->copyState(mCurrentRoadmap[u].v_state->state, oldState->state);

      //Put everything back
      for(auto nbrweights : oldEdgeMap)
      {
        Vertex nbr{nbrweights.first};
        std::pair<double,double> weights{nbrweights.second};

        std::pair<Edge,bool> current_edge = boost::edge(u,nbr,mCurrentRoadmap);
        if(!current_edge.second){
          throw ompl::Exception("Edge earlier modified no longer exists! Iterator stability issue");
        }

        mCurrentRoadmap[current_edge.first].distance = weights.first;
        mCurrentRoadmap[current_edge.first].collMeasure = weights.second;
      }

      // And now remove new edges
      for(Edge e : newEdgesToRemove)
      {
        boost::remove_edge(e,mCurrentRoadmap);
      }
    }

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
      impScoreSum += mVertexImportance[i];
    }

    // TODO - Sign here determines softmax or softmin
    double softmaxScoreSum{0.0};
    for(unsigned int i=2u ; i < mNumCurrVertices; i++)
    {
      // Prefer higher importance
      softmaxScores[i-2] = std::exp(mVertexImportance[i]/impScoreSum);
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
    std::cout<<"Initial Batch added"<<std::endl;
    // Compute initial score
    mCurrRoadmapScore = getParetoConvexHullScore(mVertexImportance, mVertexToAlphaRangeMap);
    std::cout<<"Current Score - "<<mCurrRoadmapScore<<std::endl;

    for(unsigned int i=0; i < mNumPerturbations; i++)
    {
      // Choose a sample based on its importance
      Vertex chosenVert{softMaxVertexSelection()};
      StateConPtr perturbedState(std::make_shared<StateCon>(mSpace));

      if(mVertexImportance[chosenVert] > 0.0){
        if(perturbVertexApproxGradImportance(chosenVert,perturbedState)){
          double improvement = implementPerturbation(chosenVert, perturbedState);
          // if(improvement > std::numeric_limits<double>::epsilon()){
          //   std::cout<<"Trial "<<i<<" : Vertex "<<chosenVert<<" improved score by "<<improvement<<std::endl;
          // }
        }
      }
      else{
        if(perturbVertexApproxGradNoImportance(chosenVert,perturbedState)){
          double improvement = implementPerturbation(chosenVert, perturbedState);
          // if(improvement > std::numeric_limits<double>::epsilon()){
          //   std::cout<<"Trial "<<i<<" : Vertex "<<chosenVert<<" of no importance improved score by "<<improvement<<std::endl;
          // }
        }
      }
    }

    std::cout<<"Final score - "<<mCurrRoadmapScore<<std::endl;
    //std::cout<<mSuccPerturbations<<" successful perturbations"<<std::endl;  
  }


private:

  // Testing - step size of perturbation
  double mPerturbSize;
  unsigned int mNumPerturbations;
  double mRandomOffsetSize;

  // DEBUGGING
  unsigned int mSuccPerturbations;

  std::mt19937 mGen;

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
  std::vector<double> mVertexImportance;
  std::map<Vertex,std::set<std::pair<double,double>>> mVertexToAlphaRangeMap;


  double mCurrRoadmapScore;


};


}
}

#endif
