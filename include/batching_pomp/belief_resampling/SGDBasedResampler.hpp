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
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/properties.hpp>
#include <boost/property_map/dynamic_property_map.hpp>
#include <ompl/base/spaces/RealVectorStateSampler.h>

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
                    double _perturbSize = 0.05,
  									double _dAlpha = 0.05,
  									double _probThreshold = 1.0)
  : BeliefInformedResampler(_planner)
  , mSpace(_planner.mSpace)
  , mFullRoadmap(_planner.full_g)
  , mCurrVertexNN(*(_planner.mVertexNN))
  , mStartVertex{_planner.getStartVertex()}
  , mGoalVertex{_planner.getGoalVertex()}
  . mPerturbSize{_perturbSize}
  , mProbThreshold{_probThreshold}
  , mSimpleFlag{false}
  {
  	unsigned int nAlphas = static_cast<unsigned int>(1.0/_dAlpha) + 1u;
  	alphaMinVector.resize(nAlphas);
  	alphaSecondMinVector.resize(nAlphas);
  	boost::tie(mCurrVertex,mLastVertex) = vertices(mFullRoadmap);
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
  	std::vector<Vertex> vertex_vector(BeliefInformedResampler::mBatchParams.first);

  	while(numVerticesAdded < BeliefInformedResampler::mBatchParams.first)
  	{
  		BeliefPoint query(mFullRoadmap[*mCurrVertex].v_state->state, 
  											BeliefInformedResampler::mSpace->getDimension(), -1.0);

  		// TODO - Add ellipse stuff
  		if(BeliefInformedResampler::mBeliefModel->estimate(query) < mProbThreshold) {
  			Vertex newVertex{boost::add_vertex(BeliefInformedResampler::BeliefInformedResampler::mCurrentRoadmap)};

  			// TODO : Print vertex descriptors and verify that they are in order

  			BeliefInformedResampler::BeliefInformedResampler::mCurrentRoadmap[newVertex].v_state = mFullRoadmap[*mCurrVertex].v_state;
  			vertex_vector[idx++] = newVertex;
  			++numVerticesAdded;
  		}
  		++mCurrVertex;
  	}

  	if(idx > 0u) {
  		vertex_vector.resize(idx);
  		mCurrVertexNN.add(vertex_vector);
  	}

  	mNumCurrVertices = boost::num_vertices(BeliefInformedResampler::BeliefInformedResampler::mCurrentRoadmap);

  	// Define edges based on mutual distance and compute edge weights
  	VertexIter vi, vi_end;
  	for(boost::tie(vi,vi_end)=vertices(BeliefInformedResampler::mCurrentRoadmap); vi != vi_end; ++vi) {
  		std::vector<Vertex> vertexNbrs;
  		mCurrVertexNN.nearestR(*vi,BeliefInformedResampler::mBatchParams.second,vertexNbrs);

  		for(Vertex nbr : vertexNbrs) {
  			if(nbr == *vi)
  				continue
  			if(!boost::edge(*vi, nbr, BeliefInformedResampler::mCurrentRoadmap)) {
  				std::pair<Edge,bool> new_edge = boost::add_edge(*vi,nbr,BeliefInformedResampler::mCurrentRoadmap);
  				BeliefInformedResampler::mCurrentRoadmap[new_edge.first].distance = mPlanner.vertexDistFun(*vi,nbr);
  				BeliefInformedResampler::mCurrentRoadmap[new_edge.first].blockedStatus = UNKNOWN;

  				// TODO : DANGEROUS GAME!!!
  				BeliefInformedResampler::mCurrentRoadmap[new_edge.first].collMeasure = mPlanner.computeEdgeCollisionMeasureNoStates(*vi, nbr);
  			}
  		}
  	}
  }


  void computeInitialScore()
  {
  	//Initialize vectors
  	M_to_come.resize(mNumCurrVertices);
    M_preds.resize(mNumCurrVertices);
    L_to_come.resize(mNumCurrVertices);
    L_preds.resize(mNumCurrVertices);
    M_to_go.resize(mNumCurrVertices);
    M_succs.resize(mNumCurrVertices);
    L_to_go.resize(mNumCurrVertices);
    L_succs.resize(mNumCurrVertices);
    vertexImportance.resize(mNumCurrVertices,0.0);

    mInfiniteCostEdges.clear();

    // RUN 4 DIJKSTRA SEARCHES WITH (Start,Goal) X (M,L)

    // First for M_to_come and L_to_come from start
    boost::dijkstra_shortest_paths(BeliefInformedResampler::mCurrentRoadmap, mStartVertex,
    															 &M_preds[0], &M_to_come[0],
    															 boost::get(EProps::collMeasure,BeliefInformedResampler::mCurrentRoadmap),
    															 boost::get(boost::vertex_index, BeliefInformedResampler::mCurrentRoadmap),
    															 std::less<double>(),
            											 boost::closed_plus<double>(std::numeric_limits<double>::max()), 
            											 std::numeric_limits<double>::max(),
            											 double(),
            											 boost::default_dijkstra_visitor());

    boost::dijkstra_shortest_paths(BeliefInformedResampler::mCurrentRoadmap, mStartVertex,
    															 &L_preds[0], &L_to_come[0],
    															 boost::get(EProps::distance,BeliefInformedResampler::mCurrentRoadmap),
    															 boost::get(boost::vertex_index, BeliefInformedResampler::mCurrentRoadmap),
    															 std::less<double>(),
            											 boost::closed_plus<double>(std::numeric_limits<double>::max()), 
            											 std::numeric_limits<double>::max(),
            											 double(),
            											 boost::default_dijkstra_visitor());

    boost::dijkstra_shortest_paths(BeliefInformedResampler::mCurrentRoadmap, mGoalVertex,
    															 &M_succs[0], &M_to_go[0],
    															 boost::get(EProps::collMeasure,BeliefInformedResampler::mCurrentRoadmap),
    															 boost::get(boost::vertex_index, BeliefInformedResampler::mCurrentRoadmap),
    															 std::less<double>(),
            											 boost::closed_plus<double>(std::numeric_limits<double>::max()), 
            											 std::numeric_limits<double>::max(),
            											 double(),
            											 boost::default_dijkstra_visitor());

    boost::dijkstra_shortest_paths(BeliefInformedResampler::mCurrentRoadmap, mGoalVertex,
    															 &L_succs[0], &L_to_go[0],
    															 boost::get(EProps::distance,BeliefInformedResampler::mCurrentRoadmap),
    															 boost::get(boost::vertex_index, BeliefInformedResampler::mCurrentRoadmap),
    															 std::less<double>(),
            											 boost::closed_plus<double>(std::numeric_limits<double>::max()), 
            											 std::numeric_limits<double>::max(),
            											 double(),
            											 boost::default_dijkstra_visitor());


    // FOR EACH ALPHA, GET VECTORS OVER VERTICES AND COMPUTE SCORES
    // ALSO COMPUTE ALPHA SECOND MIN VECTOR
    unsigned int nAlphas{alphaMinVector.size()};

    for(unsigned int i = 0u; i < nAlphas; i++) {

    	double currAlpha{static_cast<double>(i)/(nAlphas-1)};

    	std::vector<double> currAlphaScores(mNumCurrVertices);

    	for(unsigned int j=0u; j< mNumCurrVertices; j++) {
    		currAlphaScores[j] = currAlpha*(L_to_come[j] + L_to_go[j]) + (1.0 - currAlpha)*(M_to_come[j] + M_to_go[j])
    	}

    	// Use custom min and second min to populate alpha vectors
      double min = std::numeric_limits<double>::infinity();
      double secondMin = std::numeric_limits<double>::infinity();

      for(unsigned int j=0u; j< mNumCurrVertices; j++) {
        if(currAlphaScores[j] < min) {
          secondMin = min;
          min = currAlphaScores[j];
        }
        else if(currAlphaScores[j] < secondMin && currAlphaScores[j] > min) {
          secondMin = currAlphaScores[j];
        }
      }

      alphaMinVector[i] = min;
      alphaSecondMinVector[i] = secondMin;
    
      // Now iterate over vertices and assign 
      for(unsigned int j=0u; j< mNumCurrVertices; j++) {
        if(std::fabs(currAlphaScores[j] - min) < std::numeric_limits<double>::epsilon()) {
          vertexImportance[j] += secondMin - min;
        }
      }
    }

    mCurrRoadmapScore = getScoreFromAlphaVector();
  }

  /// The score is basically the sum of the
  /// vector divided by the number of unique values (paths)
  double getScoreFromAlphaVector(const std::vector<double>& _alphaMinVector) const
  {
    double currVal{_alphaMinVector[0]}
    double score{currVal};
    unsigned int numUniqueVals{1u};
    unsigned int nAlphas{_alphaMinVector.size()};

    for(unsigned int i = 1u; i < nAlphas; i++) {
      score += _alphaMinVector[i];
      if(std::fabs(_alphaMinVector[i] - currVal) > std::numeric_limits<double>::epsilon()) {
        numUniqueVals++;
        currVal = _alphaMinVector[i];
      }
    }

    score /= numUniqueVals;
    return score;
  }


  StateConPtr perturbVertex(const Vertex& u) const
  {
    StateConPtr perturbedState(std::make_shared<StateCon>());

    RealVectorStateSampler rvSampler(mSpace);

    rvSampler.sampleUniformNear(perturbedState->state, g[u].v_state->state, mStepSize);

    return perturbedState;
  }


  double implementPerturbation(Vertex& u, const StateConPtr& _perturbedState) const
  {
    // Use this as a proxy for new vertex
    Vertex tempNewVertex{boost::add_vertex(BeliefInformedResampler::mCurrentRoadmap)};
    g[tempNewVertex].v_state = _perturbedState;

    // To recompute old info
    std::map<Vertex,std::pair<double,double>> newNbrVertexWeights;

    // First iterate through current edges and make inf cost those that are invalidated
    OutEdgeIter ei, ei_end;
    for (boost::tie(ei,ei_end)=out_edges(u,BeliefInformedResampler::mCurrentRoadmap); ei!=ei_end; ++ei) {

      double edgeDist{mPlanner.vertexDistFun(boost::target(*ei,BeliefInformedResampler::mCurrentRoadmap), tempNewVertex)};
      
      // Edge invalid in temp new graph
      if(edgeDist > BeliefInformedResampler::mBatchParams.second) {  
        newNbrVertexWeights.insert(std::make_pair(boost::target(*ei,BeliefInformedResampler::mCurrentRoadmap),
          std::make_pair(std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity())));
      }
      else {
        // Store modified weights
        newNbrVertexWeights.insert(std::make_pair(boost::target(*ei,BeliefInformedResampler::mCurrentRoadmap),
          std::make_pair(edgeDist, mPlanner.computeEdgeCollisionMeasureNoStates(
            boost::target(*ei,BeliefInformedResampler::mCurrentRoadmap) , tempNewVertex))));
      }
    }

    std::vector<Vertex> tempVertexNbrs;
    mCurrVertexNN.nearestR(tempNewVertex,BeliefInformedResampler::mBatchParams.second,tempVertexNbrs);

    double min_M_tc{std::numeric_limits<double>::infinity()};
    double min_M_tg{std::numeric_limits<double>::infinity()};
    double min_L_tc{std::numeric_limits<double>::infinity()};
    double min_L_tg{std::numeric_limits<double>::infinity()};
    Vertex M_pred,M_succ,L_pred,L_succ;

    for(Vertex tvnbr : tempVertexNbrs) {

      std::pair<double,double> weights{std::make_pair(mPlanner.vertexDistFun(u,tvnbr), 
            mPlanner.computeEdgeCollisionMeasureNoStates(u,tvnbr))};

      if(!boost::edge(u,tvnbr,BeliefInformedResampler::mCurrentRoadmap).second) {
        newNbrVertexWeights.insert(std::make_pair(tvnbr,weights));
      }

      //Either way, update the M_tc,M_tg etc
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

    std::vector<double> tempAlphaVector(alphaMinVector);
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

        for(boost::tie(ei,ei_end)=out_edges(nbr); ei != ei_end; ++ei) {

          Vertex candidate{boost::target(*ei,BeliefInformedResampler::mCurrentRoadmap)};
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
            length = BeliefInformedResampler::mCurrentRoadmap[*ei].distance;
            meas = BeliefInformedResampler::mCurrentRoadmap[*ei].collMeasure;
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
          v_Min_L_tc = ltc + weights.first;
          v_L_pred = u;
        }
        if(min_M_tc + weights.second < M_to_come[nvnbr]) {
          v_Min_M_tc = mtc + weights.second;
          v_M_pred = u;
        }
        if(min_L_tg + weights.first < L_to_go[nvnbr]) {
          v_Min_L_tg = ltg + weights.first;
          v_L_succ = u;
        }
        if(min_M_tg + weights.second < M_to_go[nvnbr]) {
          v_Min_M_tg = mtg + weights.second;
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
                          v_Min_L_tg, v_L_succ v_Min_M_tg, v_M_succ)
        tempVertexInfoMap.insert(std::make_pair(nvnbr, nvnbr_tuple));
        tempVertexImportanceMap.insert(std::make_pair(nvnbr, vertexImportance[nvnbr]));
      }
    }


    // Now run through alphas and update score and importance for 
    unsigned int nAlphas{alphaMinVector.size()};

    for(unsigned int i = 0u; i < nAlphas; i++) {

      double currAlpha{static_cast<double>(i)/(nAlphas-1)};

      for(auto tvinfo : tempVertexInfoMap) {

        double ltc,mtc,ltg,mtg;
        Vertex lpred,mpred,lsucc,msucc;
        std::tie(ltc,lpred,mtc,mpred,ltg,lsucc,mtg,msucc) = tvinfo.second;

        double alphaScore{currAlpha*(ltc+ltg) + (1.0 - currAlpha)*(mtc+mtg)};

        // Subtract importance if it gets worse
        tempVertexImportanceMap[tvinfo.first] += alphaMinVector[i] - alphaScore;

        if( alphaScore < alphaMinVector[i]) {
          tempAlphaVector[i] = alphaScore;
        }
        else if (alphaScore < alphaSecondMinVector[i]) {
          tempSecondAlphaVector[i] = alphaScore;
        }
      }
    }

    double newRoadmapScore{getScoreFromAlphaVector(tempAlphaVector)};
    double improvement{mCurrRoadmapScore - newRoadmapScore};

    // Check if roadmap score improves
    if(newRoadmapScore < mCurrRoadmapScore) {

      mCurrRoadmapScore = newRoadmapScore;

      // Change state of vertex and other values
      mSpace->copyState(BeliefInformedResampler::mCurrentRoadmap[u].v_state->state, _perturbedState->state);
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
        std::pair<Edge,bool> potential_edge = boost::edge(u,nbr,BeliefInformedResampler::mCurrentRoadmap);

        if(potential_edge.second) {
          BeliefInformedResampler::mCurrentRoadmap[potential_edge].distance = weights.first;
          BeliefInformedResampler::mCurrentRoadmap[potential_edge].collMeasure = weights.second;

          if(weights.first == std::numeric_limits<double>::infinity()) {
            mInfiniteCostEdges.insert(potential_edge.first);
          }
        }
        else{
          std::pair<Edge,bool> new_edge = boost::add_edge(u,nbr,BeliefInformedResampler::mCurrentRoadmap);
          BeliefInformedResampler::mCurrentRoadmap[new_edge].distance = weights.first;
          BeliefInformedResampler::mCurrentRoadmap[new_edge].collMeasure = weights.second;
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
      alphaSecondMinVector = tempSecondAlphaVector;
    }
    
    // Either way, remove new vertex
    boost::remove_vertex(tempNewVertex, BeliefInformedResampler::mCurrentRoadmap);

    // Return numerical value of improvement
    return improvement;

  }


// TODO : Remove all edges with either weight infinite in main loop

  void updateRoadmap()
  {
    // First add initial samples


    // Compute initial score


    // For K times
    // Perturb and implement perturbation
  }


private:

  const ompl::base::StateSpacePtr mSpace;

  // Testing - step size of perturbation
  double mPerturbSize;

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
  std::vector<double> M_to_come;
  std::vector<Vertex> M_preds;
  std::vector<double> L_to_come;
  std::vector<Vertex> L_preds;
  std::vector<double> M_to_go;
  std::vector<Vertex> M_succs;
  std::vector<double> L_to_go;
  std::vector<Vertex> L_succs;
  std::vector<double> vertexImportance;
  std::vector<double> alphaMinVector;
  std::vector<double> alphaSecondMinVector;

  double mCurrRoadmapScore;


};


}
}

#endif