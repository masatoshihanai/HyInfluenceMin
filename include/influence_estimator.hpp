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
#include <math.h>
#include <list>

#include "hygraph.hpp"
#include "wtime.hpp"

#ifdef HAS_OMP
#include <omp.h>
#endif

class RstRRInfEstimator {
  uint64_t choose(uint64_t n, uint64_t k) {
    if (k == 0) return 1;
    return (n * choose(n - 1, k - 1)) / k;
  }

  const double apprxFctr = 1.0 - 1 / 2.71828;
  VertID n_;
  VertID LB_;
  int k_;
  uint64_t lambda0_;
  uint64_t lambda1_;
  uint64_t curNumSample_;
  double eps_;

  HyGraph* hygraph_;
  struct RstRRpath {
    std::deque<VertID> vertice_; // from v
    std::deque<HyEdgeID> hyedges_;
  };
  std::vector<RstRRpath> rstRRpaths;
  uint64_t numSample = 0;
  std::vector<std::vector<HyEdgeID>> hyEdge2RstRRpath;
  std::vector<uint64_t> deltaReduction;
  std::vector<double> vThresholds;

  void genRandomRstRRpath() {
    RstRRpath rstRRpath;
    VertID parent = -1;
    uint64_t dstFromV = 0;
    VertID v = xoshiro256p::next() % hygraph_->numVert();
    std::unordered_set<VertID> visited;
    visited.insert(v);
    std::unordered_map<HyEdgeID, uint64_t> visitedHyEdges; // <HyEdgeID, dst from v>
    rstRRpath.vertice_.push_back(v);

    /* Backward traversal from v */
    while (true) {
      double rand = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
      bool finish = false;
      for (auto hyEid: hygraph_->neighbors(v)) {
        for (auto u: hygraph_->getHyEdge(hyEid).vertices_) {
          if (u == v) continue;
          rand -= hygraph_->getIF(u,hyEid);
          if (rand < 0) {
            if (visited.count(u) == 0) {
              if (visitedHyEdges.count(hyEid) == 0) {
                visitedHyEdges.emplace(hyEid, dstFromV);
                hyEdge2RstRRpath.at(hyEid).push_back(numSample);
              }
              visited.insert(u);
              parent = u;
              rstRRpath.vertice_.push_back(u);
              rstRRpath.hyedges_.push_back(hyEid);
            }
            finish = true;
            break;
          }
        }
        if (finish) break;
      }

      if (parent == -1) {
        break;
      } else {
        v = parent;
        parent = -1;
      }
      ++dstFromV;
    }
    for (auto x: visitedHyEdges) {
      deltaReduction.at(x.first) += (dstFromV - x.second);
      //todo std::cout << "deltaReduction at " << x.first << " is " << deltaReduction.at(x.first) << " dstFromV " << dstFromV  << " x.second " << x.second << std::endl;
    }
    rstRRpaths.push_back(rstRRpath);
    ++numSample;

    // todo remove
//    std::unordered_set<VertID> s;
//    for (auto x: rstRRpath.vertice_) {
//      std::cout << x << " ";
//      s.insert(x);
//    }
//    std::cout << std::endl;
//
//    for (auto x: rstRRpath.hyedges_) {
//      std::cout << x << " ";
//    }
//    std::cout << std::endl;

  };

 public:
  void init(HyGraph* hygraph, int k, double l = 1, double eps = 0.1) {
    // todo
    xoshiro256p::initSeed(0);
    vThresholds.resize(hygraph->numVert());
    for (auto& x: vThresholds) { x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next()); }

    eps_ = eps;
    hygraph_ = hygraph;
    k_ = k;
    n_ = hygraph_->numVert();
    l = l + log(2)/log(n_); // todo
    uint64_t E_C_k = choose(hygraph_->numHyEdges(), k);
    lambda0_ = 2*n_*(apprxFctr * sqrt(l*log(n_) + log(2)) + apprxFctr*(log(E_C_k) + l*log(n_) + log(2)))/(eps*eps);
    lambda1_ = ((double) 1.0+sqrt(2)*eps/3)*(log(E_C_k) + l*log(n_) + log(log2(n_)))*n_/(eps*eps);
    deltaReduction.resize(hygraph_->numHyEdges());
    std::fill(deltaReduction.begin(), deltaReduction.end(), 0);
    hyEdge2RstRRpath.resize(hygraph_->numHyEdges());
  };

  VertID getNumSample(VertID pow2_i) {
    return pow2_i * lambda0_ / n_;
  };

  void addSampleUntil(VertID theta) {
    while (numSample < theta) {
      genRandomRstRRpath();
    }
  };

  VertID sumReduction = 0;
  std::unordered_set<HyEdgeID> topKHyEdges;
  bool runGreedy(uint64_t pow2_i = 1) {
    sumReduction = 0;
    topKHyEdges.clear();
    std::vector<long> deltaReductionTmp = std::vector<long>(hygraph_->numHyEdges(),0);
    std::vector<long> rstRRpathsTmpSize = std::vector<long>(rstRRpaths.size());
    uint64_t av = 0;
    for (int sample = 0; sample < rstRRpaths.size(); ++sample) {
      rstRRpathsTmpSize.at(sample) = rstRRpaths.at(sample).hyedges_.size();
      av += rstRRpaths.at(sample).vertice_.size();
    }
    std::cout << "average path size " << (double) av/ rstRRpaths.size() << std::endl;

    for (int i = 0; i < k_; ++i) {
      uint64_t maxReduction = 0; HyEdgeID maxEdge;
      for (HyEdgeID hyEdgeId = 0; hyEdgeId < hygraph_->numHyEdges(); ++hyEdgeId) {
        if (topKHyEdges.count(hyEdgeId) > 0) continue;
        if (maxReduction < deltaReduction.at(hyEdgeId) + deltaReductionTmp.at(hyEdgeId)) {
          if (deltaReduction.at(hyEdgeId) < (-1)*deltaReductionTmp.at(hyEdgeId)) {
            std::cout << hyEdgeId << "deltaReduction.at(hyEdgeId)" << deltaReduction.at(hyEdgeId) << "deltaReductionTmp.at(hyEdgeId)" << deltaReductionTmp.at(hyEdgeId) << std::endl;
          }
          maxReduction = deltaReduction.at(hyEdgeId) + deltaReductionTmp.at(hyEdgeId);
          maxEdge = hyEdgeId;
        }
      }
      std::cout << "num sample " << numSample << " " << i << "-th Local Max is " << maxEdge << " " << maxReduction << std::endl;
      topKHyEdges.insert(maxEdge);
      /* Update RstRR path*/
      for (uint64_t ii = 0; ii < hyEdge2RstRRpath.at(maxEdge).size(); ++ii) {
        uint64_t sampleID = hyEdge2RstRRpath.at(maxEdge).at(ii);
        RstRRpath& rstRRpath = rstRRpaths.at(sampleID);
        std::unordered_set<HyEdgeID> updated;
        bool updating = false;
        uint64_t newSize = 0;
        for (uint64_t j = 0; j < rstRRpathsTmpSize.at(sampleID); ++j) {
          HyEdgeID& hyEdgeID = rstRRpath.hyedges_.at(j);
          if (hyEdgeID == maxEdge && updating == false) {
            updating = true;
            newSize = j;
          }
          if (updating && updated.count(hyEdgeID) == 0) {
//            if (hyEdgeID == 32) {
//              std::cout << j << " " << sampleID <<  " 32 " << deltaReductionTmp.at(hyEdgeID) << " " << rstRRpathsTmpSize.at(sampleID) << " " << deltaReduction.at(hyEdgeID) << " " << rstRRpathsTmpSize.at(sampleID) - j << std::endl;
//            }
            deltaReductionTmp.at(hyEdgeID) -= rstRRpathsTmpSize.at(sampleID) - j;
            updated.insert(hyEdgeID);
          }
        }
        rstRRpathsTmpSize.at(sampleID) = newSize;
      }
    }

// todo fix
    for (auto x: topKHyEdges) {
      sumReduction += deltaReduction.at(x) + deltaReductionTmp.at(x);
    }

    std::cout << "sumReduction" << sumReduction << " per sample " << (double) sumReduction/numSample << std::endl;
    std::cout << "eps" << eps_ << " pow2_i" << pow2_i << std::endl;
    return (double) sumReduction/numSample > (double)(1.0+sqrt(2)*eps_)/pow2_i;
  };

  void evaluate(std::unordered_set<VertID>& rstVertices) {
    sumReduction = 0;
    topKHyEdges.clear();
    std::vector<long> deltaReductionTmp = std::vector<long>(hygraph_->numHyEdges(),0);
    std::vector<long> rstRRpathsTmpSize = std::vector<long>(rstRRpaths.size());
    uint64_t av = 0;
    for (int sample = 0; sample < rstRRpaths.size(); ++sample) {
      rstRRpathsTmpSize.at(sample) = rstRRpaths.at(sample).hyedges_.size();
      av += rstRRpaths.at(sample).vertice_.size();
    }
    std::cout << "average path size " << (double) av/ rstRRpaths.size() << std::endl;

    for (int i = 0; i < k_; ++i) {
      uint64_t maxReduction = 0; HyEdgeID maxEdge;
//      for (HyEdgeID hyEdgeId = 0; hyEdgeId < hygraph_->numHyEdges(); ++hyEdgeId) {
//        if (topKHyEdges.count(hyEdgeId) > 0) continue;
//        if (maxReduction < deltaReduction.at(hyEdgeId) + deltaReductionTmp.at(hyEdgeId)) {
//          if (deltaReduction.at(hyEdgeId) < (-1)*deltaReductionTmp.at(hyEdgeId)) {
//            std::cout << hyEdgeId << "deltaReduction.at(hyEdgeId)" << deltaReduction.at(hyEdgeId) << "deltaReductionTmp.at(hyEdgeId)" << deltaReductionTmp.at(hyEdgeId) << std::endl;
//          }
//          maxReduction = deltaReduction.at(hyEdgeId) + deltaReductionTmp.at(hyEdgeId);
//          maxEdge = hyEdgeId;
//        }
//      }
//      std::cout << "num sample " << numSample << " " << i << "-th Local Max is " << maxEdge << " " << maxReduction << std::endl;
      auto it = rstVertices.begin();
      maxEdge = *it;
      rstVertices.erase(it);
      topKHyEdges.insert(maxEdge);
      /* Update RstRR path*/
      for (uint64_t ii = 0; ii < hyEdge2RstRRpath.at(maxEdge).size(); ++ii) {
        uint64_t sampleID = hyEdge2RstRRpath.at(maxEdge).at(ii);
        RstRRpath& rstRRpath = rstRRpaths.at(sampleID);
        std::unordered_set<HyEdgeID> updated;
        bool updating = false;
        uint64_t newSize = 0;
        for (uint64_t j = 0; j < rstRRpathsTmpSize.at(sampleID); ++j) {
          HyEdgeID& hyEdgeID = rstRRpath.hyedges_.at(j);
          if (hyEdgeID == maxEdge && updating == false) {
            updating = true;
            newSize = j;
          }
          if (updating && updated.count(hyEdgeID) == 0) {
//            if (hyEdgeID == 32) {
//              std::cout << j << " " << sampleID <<  " 32 " << deltaReductionTmp.at(hyEdgeID) << " " << rstRRpathsTmpSize.at(sampleID) << " " << deltaReduction.at(hyEdgeID) << " " << rstRRpathsTmpSize.at(sampleID) - j << std::endl;
//            }
            deltaReductionTmp.at(hyEdgeID) -= rstRRpathsTmpSize.at(sampleID) - j;
            updated.insert(hyEdgeID);
          }
        }
        rstRRpathsTmpSize.at(sampleID) = newSize;
      }
    }

// todo fix
    for (auto x: topKHyEdges) {
      sumReduction += deltaReduction.at(x) + deltaReductionTmp.at(x);
    }
    std::cout << "sumReduction" << sumReduction << " per sample " << (double) sumReduction/numSample << std::endl;
  }

  VertID getSamplingSize () {
    return sumReduction/numSample*hygraph_->numVert()/((double) 1+sqrt(2)*eps_);
  };
};

class LiveEdgeBasedInfEstimator {
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
  HyGraph* hygraph_;
  VertID numInfluence_;
  int numSample_;

 public:
  void init(HyGraph* hyGraph, int numSample = 1) {
    xoshiro256p::initSeed(1);
    /* Construct live-edge tree */
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
    for (int i = 0; i < liveEdgeTrees_.size(); ++i) {
      for (auto& x: liveEdgeTrees_.at(i)->r_) {
        numInfluence_ += x;
      }
    }
  };

  VertID numInf() {
    return numInfluence_;
  }

  VertID estimateNumInf(HyEdgeID removeHyEdge) {
    VertID totalReduce = 0;
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
              dstVertices.push_back(localU);
            }
          }
        }
      }

      std::sort(dstVertices.begin(), dstVertices.end());

      uint64_t reducedVertices = 0;
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

#ifdef HAS_OMP
#pragma omp critical
#endif
      {
        totalReduce += reducedVertices;
      }
    }
    return numInfluence_ - totalReduce;//(numInfluence_ - totalReduce)/numSample_/hygraph_->numVert();
  };

  void update(HyEdgeID removalHyEdge) {
    liveEdgeTrees_.clear();
    hygraph_->restrictHyEdge(removalHyEdge);
    init(hygraph_,10);
  };

 private:
  void constructLiveEdgeTrees(HyGraph* hyGraph, std::vector<LiveEdgeTree*>& let) {
    std::vector<std::vector<VertIDHyEdgeID>> outedges(hyGraph->numVert());
    std::vector<bool> roots(hyGraph->numVert(), true);

    /* Construct live-edge graph */
    for (VertID v = 0; v < hyGraph->numVert(); ++v) {
      VertID parent = -1;
      bool isRoot = true;
      double rand = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
      bool finish = false;
      for (auto hyEid: hyGraph->neighbors(v)) {
        for (auto u: hygraph_->getHyEdge(hyEid).vertices_) {
          if (u == v) continue;
          if (rand < hygraph_->getIF(u,hyEid)) {
            parent = u;
            roots.at(v) = false;
            outedges.at(parent).emplace_back(v,hyEid);
            finish = true;
            break;
          } else {
            rand -= hygraph_->getIF(u,hyEid);
          }
        }
        if (finish) break;
      }
    }

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
  }
};

class MonteCarloSimInfEstimator {
  HyGraph* hyGraph_;
  std::vector<double> vThresholds;
 public:
  void init(HyGraph* hyGraph) {
    hyGraph_ = hyGraph;
    xoshiro256p::initSeed(0);
    vThresholds.resize(hyGraph_->numVert());
    for (auto& x: vThresholds) { x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next()); }
  };

  VertID estimateNumInf(HyEdgeID cuttingHedge, int numSampling = 10) {
    std::unordered_set<HyEdgeID> cut;
    cut.insert(cuttingHedge);
    return estimateNumInf(cut, numSampling);
  };

  VertID estimateNumInf(const std::unordered_set<HyEdgeID>& cuttingHedges,
                        int numSampling = 100) {
    xoshiro256p::initSeed(0);
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
      // todo
      for (auto& x: vThresholds) { x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next()); }
      std::vector<bool> infectedV(hyGraph_->numVert(), false);
      std::vector<bool> visitV(hyGraph_->numVert(), false);
      std::vector<double> vThresholdsVar(hyGraph_->numVert());
//      vThresholdsVar = vThresholds;
      double prepTime = 0;
//      for (VertID j = 0; j < hyGraph_->numVert()/1000; ++j){
//        VertID seed = xoshiro256p::next() % hyGraph_->numVert();
      for (VertID j = 0; j < hyGraph_->numVert(); ++j){
        VertID seed = j;
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
              vThresholdsVar.at(neV) -= hyGraph_->getIF(curr, neHyEdgeID);
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
    std::cout << "sum infected " << sumInfected << std::endl;
    return sumInfected / numSampling;
  };
};

#endif //HYMINSOLVER_INFLUENCE_ESTIMATOR_HPP
