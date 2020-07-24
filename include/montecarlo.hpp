/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_MONTECARLO_HPP
#define HYMINSOLVER_MONTECARLO_HPP

#include "hygraph.hpp"
#include "wtime.hpp"

class MonteCarloSimInfEstimator {
 public:
  VertID estimateNumInf(HyGraph* hyGraph, const std::unordered_set<HyEdgeID>& cuttingHedges,
                        int numSampling = 3000) {
    //xoshiro256p::initSeed(1);
    std::cout << "  Start Monte Carlo Simulation. # Sampling: " << numSampling << " ... " << std::endl;
    std::vector<VertID> numInfected(numSampling, 0);

#ifdef HAS_OMP
    int numThreads = omp_get_num_threads();
    std::vector<std::vector<double>> vertThrInitVec(numThreads);
    std::vector<std::vector<double>> vertThrVec(numThreads);

    #pragma omp parallel for
    for (int i = 0; i < numThreads; ++i) {
      vertThrInitVec.at(i) = std::vector<double>(hyGraph->numVert());
      vertThrVec.at(i) = std::vector<double>(hyGraph->numVert());
    }

    #pragma omp parallel for
    for (int i = 0; i < numSampling; ++i) {
      int thrNum = omp_get_thread_num();
      for (int j = 0; j < thrNum; ++j) { xoshiro256p::jump(); }
      std::vector<double>& vertThrInit = vertThrInitVec.at(thrNum);
      std::vector<double>& vertThr = vertThrVec.at(thrNum);
#else
    std::vector<double> vertThrInit(hyGraph->numVert());
    std::vector<double> vertThr(hyGraph->numVert());
    for (int i = 0; i < numSampling; ++i) {
#endif
      double start = getTime();
      std::vector<bool> infectedV(hyGraph->numVert(), false);
      std::vector<bool> visitV(hyGraph->numVert(), false);
      for (auto& x: vertThrInit) { x = hyGraph->genRandThr(); }

      /* LT spread from seed */
      for (VertID seed = 0; seed < hyGraph->numVert(); ++seed){
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
          for (auto& neHyEdgeID: hyGraph->neighbors(curr)) {
            HyEdge& hyEdge = hyGraph->getHyEdge(neHyEdgeID);
            for (auto& neV: hyEdge.vertices_) {
              if (cuttingHedges.count(neHyEdgeID) > 0 && hyGraph->isRestricted(neHyEdgeID, curr, neV)) continue;
              if (neV == curr) continue;
              if (infectedV.at(neV)) continue;
              if (!visitV.at(neV)) {
                vertThr.at(neV) = vertThrInit.at(neV);
                visitV.at(neV) = true;
              }
              if (vertThr.at(neV) < hyGraph->getIF(curr, neV, neHyEdgeID)) {
                infectedV.at(neV) = true;
                ++numInfected.at(i);
                frontiers.push_back(neV);
              } else {
                vertThr.at(neV) -= hyGraph->getIF(curr, neV, neHyEdgeID);
              }
            }
          }
        }
      }
      double end = getTime();
      VertID sumInfected = 0;
      VertID zeroInfected = 0;
      for (auto x: numInfected) {
        sumInfected += x;
      }
      std::cout << "\r    " << i << "-th of " << numSampling << "-th Run Finish. # infection " << numInfected.at(i) << " " << (double)(end - start) << " (sec). Avg. " << sumInfected / (i+1) << "    " <<std::flush;
    }
    std::cout << std::endl;
    VertID sumInfected = 0;
    for (auto x: numInfected) { sumInfected += x; }
    std::cout << "  Finish to Run " << numSampling << " Simulations. Total # Inf is " << sumInfected << ". Avg is " << sumInfected / numSampling << std::endl;
    return sumInfected / numSampling;
  };
};

#endif //HYMINSOLVER_MONTECARLO_HPP
