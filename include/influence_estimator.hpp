/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_INFLUENCE_ESTIMATOR_HPP
#define HYMINSOLVER_INFLUENCE_ESTIMATOR_HPP

#include <algorithm>
#include <unordered_map>

#include "hygraph.hpp"
#include "wtime.hpp"

#ifdef HAS_OMP
#include <omp.h>
#endif

class IncrementalInfEstimator {
  using VertIDHyEdgeID = std::pair<VertID,HyEdgeID>;
  using EdgeID = uint64_t;
  struct Edge {
    VertID src_;
    VertID dst_;
  };

  struct LiveEdgeTree {
    VertID root_;
    std::unordered_map<VertID, int> globalID2localID_;
    std::vector<std::vector<VertIDHyEdgeID>> outVertices_;
    std::vector<VertIDHyEdgeID> inVertices_;
    std::vector<VertID> r_;
  };

  std::vector<LiveEdgeTree*> liveEdgeTrees_;

  void constructLiveEdgeTrees(HyGraph* hyGraph, std::vector<LiveEdgeTree*>& let) {
    std::vector<std::vector<VertIDHyEdgeID>> outedges(hyGraph->numVert());
    std::vector<bool> roots(hyGraph->numVert(), true);

    std::vector<VertID> in82vert;// todo
    /* Construct live-edge graph */
    for (VertID v = 0; v < hyGraph->numVert(); ++v) {
      VertID parent = -1;
      HyEdgeID selectedHyEID = -1;
      bool isRoot = true;
      double rand = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
      for (auto hyEid: hyGraph->neighbors(v)) {
        HyEdge &hyE = hyGraph->getHyEdge(hyEid);
        if (rand < hyE.IF_) {
          selectedHyEID = hyEid;
          // todo remove
//          std::cout << "selected hyedge " << selectedHyEID << std::endl;
          break;
        } else {
          rand -= hyE.IF_;
        }
      }
      rand = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
      if (selectedHyEID != -1) {
        for (auto u: hyGraph->getHyEdge(selectedHyEID).vertices_) {
          if (u == v) continue;
          if (rand < hyGraph->getVertIF(u)) {
            parent = u;
            roots.at(v) = false;
            if (selectedHyEID == 82) {
              //std::cout << "select " << parent << " to " << v << std::endl;
              in82vert.push_back(v);
            } //todo remove
            outedges.at(parent).emplace_back(v,selectedHyEID);
            break;
          } else {
            rand -= hyGraph->getVertIF(u);
          }
        }
      }
    }
//    for (int i = 0; i < roots.size(); ++i) {
//      if (roots.at(i) == true) std::cout << i << " is root" << std::endl;
//    }
//    std::cout << "outedges" << outedges.at(0).size() << std::endl;

    /* Decompose to multiple live-edge trees */
    for (VertID v = 0; v < roots.size(); ++v) {
      if (roots.at(v) && outedges.at(v).size() != 0) {
        int localID = 0;
        LiveEdgeTree* l = new LiveEdgeTree();
        l->root_ = localID;
        l->globalID2localID_.emplace(v, localID);
        ++localID;
        l->inVertices_.emplace_back(-1,-1); // no parent from the root
        l->outVertices_.push_back(std::vector<VertIDHyEdgeID>());
        /* BFS traversal */
        std::deque<VertID> queue;
        queue.push_back(v);
        for (uint64_t i = 0; i < queue.size(); ++i) {
          VertID cur = queue.at(i);
          for (auto ne: outedges.at(cur)) {
            l->globalID2localID_.emplace(ne.first, localID);
            ++localID;
            if (l->globalID2localID_.at(ne.first) <= l->inVertices_.size()) {
              l->inVertices_.emplace_back(l->globalID2localID_.at(cur), ne.second);
              l->outVertices_.push_back(std::vector<VertIDHyEdgeID>());
            }
            l->outVertices_.at(l->globalID2localID_.at(cur)).emplace_back(l->globalID2localID_.at(ne.first), ne.second);
            queue.push_back(ne.first);
          }
        }
        /* Backwards traversal */
        l->r_ = std::vector<VertID>(localID, 0);
        while (!queue.empty()) {
          VertID cur = queue.back();
          queue.pop_back();
          for (auto ne: outedges.at(cur)) {
            l->r_.at(l->globalID2localID_.at(cur)) += l->r_.at(l->globalID2localID_.at(ne.first)) + 1;
          }
        }
        let.push_back(l);
      }
    }
    // todo remove
//    for (auto l: let) {
//      for (auto x: in82vert) {
//        if (l->globalID2localID_.count(x) > 0) {
//          std::cout << x << " r value " << l->r_.at(l->globalID2localID_.at(x)) << std::endl;
//        }
//      }
//    }
  }
  HyGraph* hygraph_;
  VertID numInfluence_;
  int numSample_;

 public:
  void init(HyGraph* hyGraph, int numSample = 1) {
    /* Construct live-edge tree */
    xoshiro256p::initSeed(0);
    hygraph_ = hyGraph;
    numSample_ = numSample;
#ifdef HAS_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < numSample; ++i) {
#ifdef HAS_OMP
      int numJump = omp_get_thread_num();
      for (int j = 0; j < numJump; ++j) { xoshiro256p::jump(); }
#endif
      std::vector<LiveEdgeTree*> newLET;
      constructLiveEdgeTrees(hyGraph, newLET);
#ifdef HAS_OMP
      #pragma omp critical
#endif
      {
        for (auto& let: newLET) {
          liveEdgeTrees_.push_back(let);
        }
      };
    }
    /* compute numInfluence */
    numInfluence_ = 0;
    std::cout << "liveEdgeTrees_.size() " << liveEdgeTrees_.size() << std::endl;
    for (int i = 0; i < liveEdgeTrees_.size(); ++i) {
      for (auto& x: liveEdgeTrees_.at(i)->r_) {
        numInfluence_ += x;
      }
    }
    std::cout << "numInf becomes: " << numInfluence_ << std::endl;
  };

  VertID estimateNumInf(HyEdgeID removeHyEdge) {
    VertID totalReduce = 0;
//    std::cout << "estimate remove " << removeHyEdge << std::endl;
#ifdef HAS_OMP
#pragma omp parallel for
#endif
    for (int i = 0; i < liveEdgeTrees_.size(); ++i) {
      LiveEdgeTree* let = liveEdgeTrees_.at(i);
      std::vector<int> dstVertices;
      for (auto& v: hygraph_->getHyEdge(removeHyEdge).vertices_) {
        if (let->globalID2localID_.count(v) == 0) continue;
        for (auto& u: hygraph_->getHyEdge(removeHyEdge).vertices_) {
          if (v == u) continue;
          if (let->globalID2localID_.count(u) == 0) continue;
          /* v and u are in the tree */
          int localV = let->globalID2localID_.at(v);
          int localU = let->globalID2localID_.at(u);
          if (localV < localU) {
            if (let->inVertices_.at(localU).first == localV
                && let->inVertices_.at(localU).second == removeHyEdge) {
              /* v to u exists */
              //std::cout << i << " remove " << v << " to " << u << std::endl;
              dstVertices.push_back(localU);
            }
          }
        }
      }

      std::sort(dstVertices.begin(), dstVertices.end());

      uint64_t reducedVertices = 0;
//      if (i == 27) { // todo
//        std::cout << "debug " << std::endl;
//      }
      while (!dstVertices.empty()) {
        VertID dst = dstVertices.back();
        dstVertices.pop_back();
        int cur = let->inVertices_.at(dst).first;
        reducedVertices += let->r_.at(dst) + 1;
        while (cur != let->root_) {
          if (std::binary_search(dstVertices.begin(), dstVertices.end(), cur)) {
            break;
          }
          cur = let->inVertices_.at(cur).first;
          reducedVertices += let->r_.at(dst) + 1;
        }
      }
//      if (i == 27) std::cout << "reduce vertices " << reducedVertices << std::endl;// todo

#ifdef HAS_OMP
#pragma omp critical
#endif
      {
        totalReduce += reducedVertices;
      }
    }
//    std::cout << "numInfluence_" << numInfluence_ << std::endl;
//    std::cout << "total_" << totalReduce << std::endl;
//    std::cout << numInfluence_ - totalReduce << std::endl;
    return numInfluence_ - totalReduce;//(numInfluence_ - totalReduce)/numSample_/hygraph_->numVert();
  };

  VertID numInf() {
    return numInfluence_;
  }

  void update(HyEdgeID removalHyEdge) {
    liveEdgeTrees_.clear();
    hygraph_->restrictHyEdge(removalHyEdge);
    init(hygraph_,10);
  };
};

class MonteCarloSimInfEstimator {
  HyGraph* hyGraph_;
 public:
  void init(HyGraph* hyGraph) {
    hyGraph_ = hyGraph;
    xoshiro256p::initSeed(0);
  };

  VertID estimateNumInf(HyEdgeID cuttingHedge, int numSampling = 10) {
    std::unordered_set<HyEdgeID> cut;
    cut.insert(cuttingHedge);
    estimateNumInf(cut, numSampling);
  };

  VertID estimateNumInf(const std::unordered_set<HyEdgeID>& cuttingHedges,
                        int numSampling = 1) {
    std::cout << "    Start Cascade Simulation. # Sampling: " << numSampling << " ... " << std::flush;
    std::vector<VertID> numInfected(numSampling, 0);
#ifdef HAS_OMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < numSampling; ++i) {
#ifdef HAS_OMP
      int numJump = omp_get_thread_num();
      for (int j = 0; j < numJump; ++j) { xoshiro256p::jump(); }
#endif
      double start = getTime();
      std::vector<bool> infectedV(hyGraph_->numVert(), false);
      std::vector<bool> visitV(hyGraph_->numVert(), false);
      std::vector<double> vThresholds(hyGraph_->numVert());
      for (auto& x: vThresholds) { x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next()); }
      std::vector<double> vThresholdsVar(hyGraph_->numVert());
      vThresholdsVar = vThresholds;
      double prepTime = 0;
      for (VertID seed = 0; seed < hyGraph_->numVert(); ++seed){
        std::fill(infectedV.begin(), infectedV.end(), false);
        std::fill(visitV.begin(), visitV.end(), false);

        /* Run HyLT Infection model */
        std::deque<VertID> frontiers;
        frontiers.push_back(seed);
        infectedV.at(seed) = true;
        visitV.at(seed) = true;
        while (frontiers.size() != 0) {
          VertID curr = frontiers.front();
          frontiers.pop_front();
          for (auto& neHyEdgeID: hyGraph_->neighbors(curr)) {
            if (cuttingHedges.count(neHyEdgeID) > 0) { continue; }
            HyEdge& hyEdge = hyGraph_->getHyEdge(neHyEdgeID);
            for (auto& neV: hyEdge.vertices_) {
              if (infectedV.at(neV)) continue;
              if (!visitV.at(neV)) {
                vThresholdsVar.at(neV) = vThresholds.at(neV);
                visitV.at(neV) = true;
              }
              vThresholdsVar.at(neV) -= hyEdge.IF_*hyGraph_->getVertIF(neV);
              if (vThresholdsVar.at(neV) < 0) {
                infectedV.at(neV) = true;
                ++numInfected.at(i);
                frontiers.push_back(neV);
              }
            }
          }
        }
      }
      double end = getTime();
      if (i == 0) std::cout << "Initial Run Finish. Expected Total Time is " << (double)(end - start)*numSampling << " (sec)" << std::flush;
    }
    std::cout << " Finish to Run " << numSampling << " Simulations." << std::endl;
    VertID sumInfected = 0;
    for (auto x: numInfected) { sumInfected += x; }
    return sumInfected / numSampling;
  };
};

#endif //HYMINSOLVER_INFLUENCE_ESTIMATOR_HPP
