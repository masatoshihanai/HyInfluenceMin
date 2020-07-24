/*
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

bool VERBOSE = false;

#include <algorithm>
#include <iostream>
#include <unistd.h>
#include <vector>
#include <utility>
#include <unordered_set>
#include <random>

#include "external/xoshiro256plus.hpp"
#include "hygraph.hpp"
#include "rstrrpath_estimator.hpp"
#include "livepath_estimator.hpp"
#include "montecarlo.hpp"
#include "wtime.hpp"

int main(int argc, char** argv) {
  xoshiro256p::initSeed(1);
  int opt;
  std::string measureType = "delete";
  bool ANALYZE = false;
  std::string repeatInterval = "daily";
  while ((opt = getopt(argc, argv, "hvm:s:ar:u:")) != -1) {
    switch (opt) {
      case 'h': {
        std::cerr << " HyMin - Hypergraph-based Influence Minimization Solver\n"
                  << "   Usage: $ " << argv[0] << " [-h] [-v] [-a] [-m <safety measure>] [-s <random seed>] [-r <repeat interval>] [-u time unit] <checkin-file> <# of restriction activities> \n"
                  << "      -h: Help\n"
                  << "      -v: Verbose. Output details\n"
                  << "      -m: Safety Measure type. delete, shrink, or split (default: delete)\n"
                  << "      -a: Analyze with MonteCarlo Simulation (default: off) !! Time-Consuming !!\n"
                  << "      -s: Random seed (default: 1)\n"
                  << "      -r: Repeat interval. weekly or daily (default: daily) \n"
                  << "      -u: Time unit for activity (hour).  (default: 12). Check-ins within x hours are considered to be the same activity. \n"
                  << "\n"
                  << "    Example to use:\n"
                  << "      $ " << argv[0] << " -v -m delete ../data/brightkite.checkin 100"
                  << std::endl;
        return 0;
      }
      case 'v': {
        VERBOSE = true;
        break;
      }
      case 'm': {
        measureType = std::string(optarg);
        break;
      }
      case 's': {
        xoshiro256p::initSeed(std::atol(optarg));
        break;
      }
      case 'a': {
        ANALYZE = true;
        break;
      }
      case 'r': {
        repeatInterval = std::string(optarg);
        break;
      }
      case 'u': {
        TIME_UNIT_HOUR = std::atoi(optarg);
        break;
      }
      default:
        return 1;
    }
  }

  if (argc - optind < 2) {
    std::cerr << "Run \"$" << argv[0] << " -h\" for help." << std::endl;
    return 1;
  }

  /* Init */
  HyGraph hygraph;
  hygraph.init(argv[optind], measureType, repeatInterval);

  int k = std::atoi(argv[optind+1]);

  std::cout << "--------- Initial Info --------" << std::endl;
  std::cout << " Data: " << argv[optind] << std::endl;
  std::cout << " # check-in:      " << hygraph.numCheckIn() << std::endl;
  std::cout << " # vertices:      " << hygraph.numVert() << std::endl;
  std::cout << " # hyperedges:    " << hygraph.numHyEdges() << std::endl;
  std::cout << " # restriction k: " << k << std::endl;
//  std::cout << " Algorithm:       " << algType << std::endl;
  std::cout << " Safety Measure:  " << measureType << std::endl;
  std::cout << " Repeat interval: " << repeatInterval << std::endl;
  std::cout << " Time unit:       " << TIME_UNIT_HOUR << " (hour)" << std::endl;
  std::cout << "-------------------------------" << std::endl;

  /* Start */
  std::cout << "==== Start Computation ====" << std::endl;
  /* RstRR-path-based */
  double start = getTime();
  HyIDSet topKHyEdges;
  RstRRPathEstimator rstRRPathEstimator;
  uint64_t theta;
  rstRRPathEstimator.init(&hygraph, k);
  for (uint64_t pow2_i = 1; pow2_i < hygraph.numVert(); pow2_i *= 2) {
    theta = rstRRPathEstimator.getNumSample(pow2_i);
    rstRRPathEstimator.addSampleUntil(theta);
    if (rstRRPathEstimator.runGreedy(topKHyEdges, pow2_i)) {
      break;
    }
  }

  theta = rstRRPathEstimator.getSamplingSize();
  rstRRPathEstimator.addSampleUntil(theta);
  rstRRPathEstimator.runGreedy(topKHyEdges);
  std::cout << "==== Finish ====" << std::endl;
  double end = getTime();
  std::cout << "Execution Time(sec): " << (double)(end - start) << std::endl;

  std::cout << "---- Start Analyze Result ----" << std::endl;
  std::cout << "  Restriction Activities (HyEdgeID): ";
  for (auto& x: topKHyEdges) { std::cout << x << ","; }
  std::cout << std::endl;

  /* Random */
  VertSet topKrandom;
  for (int i = 0; i < k; ++i) {
    HyEdgeID e;
    while (topKrandom.count((e = (xoshiro256p::next() % hygraph.numHyEdges()))) != 0) {};
    topKrandom.insert(e);
  }
  std::cout << "  Restriction Activities For Random: ";
  for (auto& x: topKrandom) { std::cout << x << ","; }
  std::cout << std::endl;

  /* Degree-based */
  VertSet topKdegree;
  using HyDeg = std::pair<HyEdgeID, double>;
  std::vector<HyDeg> degrees;
  for (HyEdgeID i = 0; i < hygraph.numHyEdges(); ++i) { degrees.emplace_back(i, hygraph.getHyEdge(i).vertices_.size()); }
  std::sort(degrees.begin(), degrees.end(), [](const HyDeg& left, const HyDeg& right){ return left.second > right.second; });
  for (int i = 0; i < k; ++i) { topKdegree.insert(degrees.at(i).first); }
  std::cout << "  Restriction Activities For Degree: ";
  for (auto& x: topKdegree) { std::cout << x << ","; }
  std::cout << std::endl;

  /* Evaluate vs random and degree */
  double numRedRandom = rstRRPathEstimator.numInfByRestriction(topKrandom);
  double numRedDegree = rstRRPathEstimator.numInfByRestriction(topKdegree);
  double numRed = rstRRPathEstimator.numInfByRestriction(topKHyEdges);
  double numInf = rstRRPathEstimator.initialNumInf();
  std::cout << "  # Inf:              " << numInf << std::endl;
  std::cout << "  # Inf in RstRRPath: " << numRed << std::endl;
  std::cout << "  # Inf in degree:    " << numRedDegree << std::endl;
  std::cout << "  # Inf in random:    " << numRedRandom << std::endl;
  std::cout << "  (# Inf in RstRRPath) / (# Inf in degree): " << numRed/numRedDegree << std::endl;
  std::cout << "  (# Inf in RstRRPath) / (# Inf in random): " << numRed/numRedRandom << std::endl;
  std::cout << "  (# Inf in RstRRPath) / (# Inf):           " << numRed/numInf << std::endl;
  std::cout << "---- Finish ----" << std::endl;

  if (ANALYZE) {
    /*--------------------------------------------*/
    /* Analyze Result with Monte Carlo Simulation */
    /*--------------------------------------------*/
    std::cout << "---- Start Analyze Result with Monte Carlo Simulation ----" << std::endl;
    /* calculate # of influence (infected vertices) */
    MonteCarloSimInfEstimator monteEstimator;
    double numInf = monteEstimator.estimateNumInf(&hygraph, topKHyEdges, 3000);
    std::cout << "  Average # Influence from all vertices: " << numInf << std::endl;

    /* Compare vs Degree-based */
    std::cout << "  << Comparing to Degree-based >>" << std::endl;

    /* Degree-based */
    using HyDeg = std::pair<HyEdgeID, uint64_t>;
    std::vector<HyDeg> degrees;
    for (HyEdgeID i = 0; i < hygraph.numHyEdges(); ++i) {
      degrees.emplace_back(i, hygraph.getHyEdge(i).vertices_.size());
    }
    std::sort(degrees.begin(), degrees.end(),
              [](const HyDeg &left, const HyDeg &right) { return left.second > right.second; });
    HyIDSet topKdegree;
    for (int i = 0; i < k; ++i) {
      topKdegree.insert(degrees.at(i).first);
    }
    MonteCarloSimInfEstimator degreeEstimator;
    double numInfDegree = degreeEstimator.estimateNumInf(&hygraph, topKdegree, 3000);
    std::cout << "  Average # Influence from all vertices: " << numInfDegree << std::endl;

    std::cout << "---- Finish ----" << std::endl;
  }

  return 0;
}