/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_LIVEPATH_ESTIMATOR_HPP
#define HYMINSOLVER_LIVEPATH_ESTIMATOR_HPP

//class LiveEdgeBasedInfEstimator {
//  using VertIDHyEdgeID = std::pair<VertID,HyEdgeID>;
//  using EdgeID = uint64_t;
//  struct Edge {
//    VertID src_;
//    VertID dst_;
//  };
//
//  struct LiveEdgeTree {
//    VertID root_;
//    std::unordered_map<VertID, int> globalID2localID_;
//    std::vector<std::vector<VertIDHyEdgeID>> outVertices_;
//    std::vector<VertIDHyEdgeID> inVertices_;
//    std::vector<VertID> r_;
//  };
//  std::vector<LiveEdgeTree*> liveEdgeTrees_;
//  HyGraph* hygraph_;
//  VertID numInfluence_;
//  int numSample_;
//  int k_ = 0;
//
// public:
//  void init(HyGraph* hyGraph, int k, int numSample) {
//    xoshiro256p::initSeed(1);
//    /* Construct live-edge tree */
//    hygraph_ = hyGraph;
//    numSample_ = numSample;
//    k_ = k;
//#ifdef HAS_OMP
//#pragma omp parallel for
//#endif
//    for (int i = 0; i < numSample; ++i) {
//#ifdef HAS_OMP
//      int numJump = omp_get_thread_num();
//      for (int j = 0; j < numJump; ++j) { xoshiro256p::jump(); }
//#endif
//      std::vector<LiveEdgeTree*> newLET;
//      constructLiveEdgeTrees(hyGraph, newLET);
//#ifdef HAS_OMP
//#pragma omp critical
//#endif
//      {
//        for (auto& let: newLET) {
//          liveEdgeTrees_.push_back(let);
//        }
//      };
//    }
//    /* compute numInfluence */
//    numInfluence_ = 0;
//    for (int i = 0; i < liveEdgeTrees_.size(); ++i) {
//      for (auto& x: liveEdgeTrees_.at(i)->r_) {
//        numInfluence_ += x;
//      }
//    }
//  };
//
//  VertID numInf() {
//    return numInfluence_;
//  }
//
//  VertID estimateNumInf(HyEdgeID removeHyEdge) {
//    VertID totalReduce = 0;
//#ifdef HAS_OMP
//#pragma omp parallel for
//#endif
//    for (int i = 0; i < liveEdgeTrees_.size(); ++i) {
//      LiveEdgeTree* let = liveEdgeTrees_.at(i);
//      std::vector<int> dstVertices;
//      for (auto& v: hygraph_->getHyEdge(removeHyEdge).vertices_) {
//        if (let->globalID2localID_.count(v) == 0) continue;
//        for (auto& u: hygraph_->getHyEdge(removeHyEdge).vertices_) {
//          if (v == u) continue;
//          if (let->globalID2localID_.count(u) == 0) continue;
//          /* v and u are in the tree */
//          int localV = let->globalID2localID_.at(v);
//          int localU = let->globalID2localID_.at(u);
//          if (localV < localU) {
//            if (let->inVertices_.at(localU).first == localV
//                && let->inVertices_.at(localU).second == removeHyEdge) {
//              /* v to u exists */
//              dstVertices.push_back(localU);
//            }
//          }
//        }
//      }
//
//      std::sort(dstVertices.begin(), dstVertices.end());
//
//      uint64_t reducedVertices = 0;
//      while (!dstVertices.empty()) {
//        VertID dst = dstVertices.back();
//        dstVertices.pop_back();
//        int cur = let->inVertices_.at(dst).first;
//        reducedVertices += let->r_.at(dst) + 1;
//        while (cur != let->root_) {
//          if (std::binary_search(dstVertices.begin(), dstVertices.end(), cur)) {
//            break;
//          }
//          cur = let->inVertices_.at(cur).first;
//          reducedVertices += let->r_.at(dst) + 1;
//        }
//      }
//
//#ifdef HAS_OMP
//#pragma omp critical
//#endif
//      {
//        totalReduce += reducedVertices;
//      }
//    }
//    return numInfluence_ - totalReduce;//(numInfluence_ - totalReduce)/numSample_/hygraph_->numVert();
//  };
//
//  void update(HyEdgeID removalHyEdge) {
//    liveEdgeTrees_.clear();
//    hygraph_->restrictHyEdge(removalHyEdge);
//    init(hygraph_,k_, numSample_);
//  };
//
//  void runGreedy(HyIDSet& topKHyEdges) {
//    for (int i = 0; i < k_; ++i) {
//      std::cout << "  Greedily Getting " << i << "-th Hyperedge to be restricted ..." << std::endl;
//      HyEdgeID e_min; VertID minInf = VID_MAX;
//      for (HyEdgeID e = 0; e < hygraph_->numHyEdges(); ++e) {
//        if (topKHyEdges.count(e) > 0) continue;
//        VertID numInf = estimateNumInf(e);
//        if (numInf < minInf) {
//          minInf = numInf;
//          e_min = e;
//        }
//      }
//      update(e_min);
//      topKHyEdges.insert(e_min);
//      std::cout << "  Finish. " << i << "-th Local Min is " << e_min << std::endl;
//    }
//  }
//
// private:
//  void constructLiveEdgeTrees(HyGraph* hyGraph, std::vector<LiveEdgeTree*>& let) {
//    std::vector<std::vector<VertIDHyEdgeID>> outedges(hyGraph->numVert());
//    std::vector<bool> roots(hyGraph->numVert(), true);
//
//    /* Construct live-edge graph */
//    for (VertID v = 0; v < hyGraph->numVert(); ++v) {
//      VertID parent = -1;
//      bool isRoot = true;
//      double rand = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
//      bool finish = false;
//      for (auto hyEid: hyGraph->neighbors(v)) {
//        for (auto u: hygraph_->getHyEdge(hyEid).vertices_) {
//          if (u == v) continue;
//          if (rand < hygraph_->getIF(u,hyEid)) {
//            parent = u;
//            roots.at(v) = false;
//            outedges.at(parent).emplace_back(v,hyEid);
//            finish = true;
//            break;
//          } else {
//            rand -= hygraph_->getIF(u,hyEid);
//          }
//        }
//        if (finish) break;
//      }
//    }
//
//    /* Decompose to multiple live-edge trees */
//    for (VertID v = 0; v < roots.size(); ++v) {
//      if (roots.at(v) && outedges.at(v).size() != 0) {
//        int localID = 0;
//        LiveEdgeTree* l = new LiveEdgeTree();
//        l->root_ = localID;
//        l->globalID2localID_.emplace(v, localID);
//        ++localID;
//        l->inVertices_.emplace_back(-1,-1); // no parent from the root
//        l->outVertices_.push_back(std::vector<VertIDHyEdgeID>());
//        /* BFS traversal */
//        std::deque<VertID> queue;
//        queue.push_back(v);
//        for (uint64_t i = 0; i < queue.size(); ++i) {
//          VertID cur = queue.at(i);
//          for (auto ne: outedges.at(cur)) {
//            l->globalID2localID_.emplace(ne.first, localID);
//            ++localID;
//            if (l->globalID2localID_.at(ne.first) <= l->inVertices_.size()) {
//              l->inVertices_.emplace_back(l->globalID2localID_.at(cur), ne.second);
//              l->outVertices_.push_back(std::vector<VertIDHyEdgeID>());
//            }
//            l->outVertices_.at(l->globalID2localID_.at(cur)).emplace_back(l->globalID2localID_.at(ne.first), ne.second);
//            queue.push_back(ne.first);
//          }
//        }
//        /* Backwards traversal */
//        l->r_ = std::vector<VertID>(localID, 0);
//        while (!queue.empty()) {
//          VertID cur = queue.back();
//          queue.pop_back();
//          for (auto ne: outedges.at(cur)) {
//            l->r_.at(l->globalID2localID_.at(cur)) += l->r_.at(l->globalID2localID_.at(ne.first)) + 1;
//          }
//        }
//        let.push_back(l);
//      }
//    }
//  }
//};

#endif //HYMINSOLVER_LIVEPATH_ESTIMATOR_HPP
