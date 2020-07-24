/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_RSTRRPATH_ESTIMATOR_HPP
#define HYMINSOLVER_RSTRRPATH_ESTIMATOR_HPP

#include <algorithm>
#include <unordered_map>
#include <math.h>
#include <list>

#include "hygraph.hpp"
#include "wtime.hpp"

#ifdef HAS_OMP
#include <omp.h>
#endif

class RstRRPathEstimator {
  double logchoose(int n, int k) {
    double ret = 0;
    for (uint64_t i = n - k + 1; i <= n; i++) { ret += log(i); }
    for (uint64_t i = 1; i <= k; i++) { ret -= log(i); }
    return ret;
  }

  /* Parameters */
  const double ONE_MINUS_1_e = 1.0 - 1.0 / 2.71828;
  VertID n_;
  VertID LB_;
  int k_;
  double lambda0_;
  double lambda1_;
  uint64_t curNumSample_;
  double eps_;
  uint64_t numSample_ = 0;

  /* Graph */
  HyGraph* hygraph_;
  struct RstRRpath {
    std::deque<VertID> vertice_; // from v
    std::deque<HyEdgeID> hyedges_;
  };
  using RstRRpaths = std::vector<RstRRpath>;

  /* Sample */
  RstRRpaths rstRRpaths_;
  using HyEID2RstRRpath = std::vector<std::vector<HyEdgeID>>;
  HyEID2RstRRpath hyEdge2RstRRpath_;
  std::vector<uint64_t> reduction_;

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
      double rand = hygraph_->genRandThr();
      bool finish = false;
      for (auto hyEid: hygraph_->neighbors(v)) {
        for (auto u: hygraph_->getHyEdge(hyEid).vertices_) {
          if (u == v) continue;
          if (rand < hygraph_->getIF(u, v, hyEid)) {
            if (visited.count(u) == 0) {
              if (visitedHyEdges.count(hyEid) == 0) {
                visitedHyEdges.emplace(hyEid, dstFromV);
                hyEdge2RstRRpath_.at(hyEid).push_back(rstRRpaths_.size());
              }
              visited.insert(u);
              parent = u;
              rstRRpath.vertice_.push_back(u);
              rstRRpath.hyedges_.push_back(hyEid);
            }
            finish = true;
            break;
          } else {
            rand -= hygraph_->getIF(u, v, hyEid);
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

    // todo fix bug
    for (auto x: visitedHyEdges) {
      reduction_.at(x.first) += (dstFromV - x.second);
    }
    if (rstRRpath.hyedges_.size() != 0) {
      rstRRpaths_.push_back(rstRRpath);
    }
  };

 public:
  void init(HyGraph* hygraph, int k, double l_0 = 1.0, double eps = 0.1) {
    //xoshiro256p::initSeed(99);
    eps_ = eps;
    hygraph_ = hygraph;
    k_ = k;
    n_ = hygraph_->numVert();
    double l = l_0 + log(2)/log(n_); // todo
    double log_E_C_k = logchoose(hygraph_->numHyEdges(), k);
    lambda0_ = 2 * n_ * (ONE_MINUS_1_e * sqrt(l * log(n_) + log(2)) + ONE_MINUS_1_e * (log_E_C_k + l * log(n_) + log(2))) / (eps * eps);
    lambda1_ = (1.0+sqrt(2)*eps/3)*(log_E_C_k + l*log(n_) + log(log2(n_)))*n_/(eps*eps);
    reduction_.resize(hygraph_->numHyEdges());
    std::fill(reduction_.begin(), reduction_.end(), 0);
    hyEdge2RstRRpath_.resize(hygraph_->numHyEdges());
    if (VERBOSE) {
      std::cout << "  ----- Parameters -----" << std::endl;
      std::cout << "    l:       " << l_0 << std::endl;
      std::cout << "    eps:     " << eps_ << std::endl;
      std::cout << "    lambda0: " << lambda0_ << std::endl;
      std::cout << "    lambda1: " << lambda1_ << std::endl;
      std::cout << "  -----------------------" << std::endl;
    }
  };

  uint64_t getNumSample(uint64_t pow2_i) {
    return pow2_i * lambda1_ / n_;
  };

  void addSampleUntil(uint64_t theta) {
    if (VERBOSE) {
      std::cout << "  Add sample from " << numSample_ << " to " << theta << std::endl;
    }
    while (numSample_ < theta) {
      genRandomRstRRpath();
      ++numSample_;
    }
    if (VERBOSE) {
      std::cout << "  Finish Sampling " << std::endl;
    }
  };
  
  VertID sumReduction = 0;
  bool runGreedy(HyIDSet& topKHyEdges, uint64_t pow2_i = 1) {
    if (VERBOSE) {
      std::cout << "  Run greedy search. # Sample: " << numSample_ << std::endl;
    }
    sumReduction = 0;
    topKHyEdges.clear();
    std::vector<long> deltaReduction = std::vector<long>(hygraph_->numHyEdges(), 0);
    std::vector<long> rstRRpathsTmpSize = std::vector<long>(rstRRpaths_.size());
    for (int sample = 0; sample < rstRRpaths_.size(); ++sample) {
      rstRRpathsTmpSize.at(sample) = rstRRpaths_.at(sample).hyedges_.size();
    }

    if (VERBOSE) {
      std::cout << "    (i-th, EdgeID, # Reduction):" << std::flush;
    }
    for (int i = 0; i < k_; ++i) {
      /* Find k-th local optimum */
      uint64_t maxReduction = 0; HyEdgeID maxEdge = 0;
      for (HyEdgeID hyEdgeId = 0; hyEdgeId < hygraph_->numHyEdges(); ++hyEdgeId) {
        if (topKHyEdges.count(hyEdgeId) > 0) continue;
        if (maxReduction <= reduction_.at(hyEdgeId) + deltaReduction.at(hyEdgeId)) {
          maxReduction = reduction_.at(hyEdgeId) + deltaReduction.at(hyEdgeId);
          maxEdge = hyEdgeId;
        }
      }
      topKHyEdges.insert(maxEdge);
      sumReduction += maxReduction;

      if (VERBOSE) {
        std::cout << "(" << i << "-th," << maxEdge << "," << maxReduction << ")" << std::flush;
      }

      /* Update RstRR path */
      for (uint64_t ii = 0; ii < hyEdge2RstRRpath_.at(maxEdge).size(); ++ii) {
        uint64_t sampleID = hyEdge2RstRRpath_.at(maxEdge).at(ii);
        RstRRpath& rstRRpath = rstRRpaths_.at(sampleID);
        /* Find index for maxEdge */
        uint64_t indxMax = 0;
        for (; indxMax < rstRRpathsTmpSize.at(sampleID); ++indxMax) {
          if (rstRRpath.hyedges_.at(indxMax) == maxEdge) break;
        }

        if (indxMax != rstRRpathsTmpSize.at(sampleID)) {
          /* Update deltaReduction */
          std::unordered_set<HyEdgeID> restricted;
          for (uint64_t j = 0; j < rstRRpathsTmpSize.at(sampleID); ++j) {
            HyEdgeID hyEdgeID = rstRRpath.hyedges_.at(j);
            VertID src = rstRRpath.vertice_.at(j);
            VertID dst = rstRRpath.vertice_.at(j+1);
            if (restricted.count(hyEdgeID) == 0) {
              if (j < indxMax) {
                if (hygraph_->isRestricted(hyEdgeID, src, dst)) {
                  deltaReduction.at(hyEdgeID) -=  rstRRpathsTmpSize.at(sampleID) - indxMax;
                  restricted.insert(hyEdgeID);
                }
              } else {
                if (hygraph_->isRestricted(hyEdgeID, src, dst)) {
                  deltaReduction.at(hyEdgeID) -= rstRRpathsTmpSize.at(sampleID) - j;
                  restricted.insert(hyEdgeID);
                }
              }
            }
          }
          rstRRpathsTmpSize.at(sampleID) = indxMax;
        }
      }
    }
    if (VERBOSE) std::cout << std::endl;

    if (VERBOSE) {
      std::cout << "  Finish running greedy search.";
      std::cout << " # reductions from all vertices: " << sumReduction*hygraph_->numVert() << std::endl;
    }
    return (double) sumReduction / numSample_ - (1.0 + sqrt(2) * eps_) / pow2_i > 0;
  };

  uint64_t getSamplingSize () {
    return lambda0_ / ((double) sumReduction / numSample_ * hygraph_->numVert()) * ((double) 1 + sqrt(2) * eps_);
  };

  double initialNumInf() {
    uint64_t ret = 0;
    uint64_t numNonZero = 0;
    for (uint64_t i = 0; i < rstRRpaths_.size(); ++i) {
      ret += rstRRpaths_.at(i).hyedges_.size();
    }
    return (double) (ret * hygraph_->numVert()) / numSample_ ;
  };

  double numInfByRestriction(const HyIDSet& topKHyEdges) {
    uint64_t sumInf = 0;
    for (int i = 0; i < rstRRpaths_.size(); ++i) {
      RstRRpath& rstRRpath = rstRRpaths_.at(i);
      uint64_t j = 0;
      for (; j < rstRRpath.hyedges_.size(); ++j) {
        HyEdgeID rstEdge = rstRRpath.hyedges_.at(j);
        VertID src = rstRRpath.vertice_.at(j);
        VertID dst = rstRRpath.vertice_.at(j+1);
        if (topKHyEdges.count(rstEdge) > 0 && hygraph_->isRestricted(rstEdge, src, dst)) {
          break;
        }
      }
      sumInf += j;
    }

    return (double) (sumInf * hygraph_->numVert()) / numSample_;
  }
};

#endif //HYMINSOLVER_RSTRRPATH_ESTIMATOR_HPP
