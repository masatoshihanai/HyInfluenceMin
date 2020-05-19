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
#include "influence_estimator.hpp"
#include "wtime.hpp"

int main(int argc, char** argv) {
  xoshiro256p::initSeed(1);

  int opt;
  std::string algType = "greedy";
  std::string measureType = "delete";
  while ((opt = getopt(argc, argv, "hva:m:s:")) != -1) {
    switch (opt) {
      case 'h': {
        std::cerr << " HyMin - Hypergraph-based Influence Minimization Solver\n"
                  << "   Usage: $ " << argv[0] << " [-h] [-v] [-t <algorithm type>] <checkin-file> <# of restriction activities> \n"
                  << "      -h: Help\n"
                  << "      -v: Verbose. Output details\n"
                  << "      -a: Algorithm type. greedy, random, or degree (default: greedy)\n"
                  << "      -m: Safety Measure type. delete, shrink, or split (default: delete)\n"
                  << "      -s: random seed (default: 1)\n"
                  << "\n"
                  << "    Example to use:\n"
                  << "      $ " << argv[0] << " -v -a greedy -m delete ../data/brightkite.checkin 100"
                  << std::endl;
        return 0;
      }
      case 'v': {
        VERBOSE = true;
        break;
      }
      case 'a': {
        algType = std::string(optarg);
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
  hygraph.init(argv[optind]);

  int k = std::atoi(argv[optind+1]);
  std::cout << "--------- Initial Info --------" << std::endl;
  std::cout << " Data: " << argv[optind] << std::endl;
  std::cout << " # vertices:      " << hygraph.numVert() << std::endl;
  std::cout << " # hyperedges:    " << hygraph.numHyEdges() << std::endl;
  std::cout << " # restriction k: " << k << std::endl;
  std::cout << " Algorithm:       " << algType << std::endl;
  std::cout << " Safety Measure:  " << measureType << std::endl;
  std::cout << "-------------------------------" << std::endl;

  /* Start */
  std::cout << "==== Start Computation ====" << std::endl;
  double start = getTime();
  std::unordered_set<HyEdgeID> topKHyEdges; // pair = <hyedge id for cutting, # of influences by this cutting>
  if (algType == "random") {
    /* random */
    for (int i = 0; i < k; ++i) {
      HyEdgeID e;
      while (topKHyEdges.count((e = (xoshiro256p::next() % hygraph.numHyEdges()))) != 0);
      topKHyEdges.insert(e);
    }
  } else if (algType == "degree") {
    IncrementalInfEstimator incEstimator;
    int numSample = 10;
    incEstimator.init(&hygraph, numSample);
    std::vector<std::pair<HyEdgeID,uint64_t>> degrees;
    for (HyEdgeID i = 0; i < hygraph.numHyEdges(); ++i) {
      degrees.emplace_back(i, hygraph.getHyEdge(i).vertices_.size());
    }
    std::sort(degrees.begin(), degrees.end(),
      [](const std::pair<HyEdgeID,uint64_t>& left, const std::pair<HyEdgeID,uint64_t>& right){ return left.second > right.second; });
    for (int i = 0; i < k; ++i) {
      std::cout << degrees.at(i).first << std::endl;
      topKHyEdges.insert(degrees.at(i).first);
      VertID numInf = incEstimator.estimateNumInf(degrees.at(i).first);
      std::cout << "numInf for " << degrees.at(i).first << " is " << numInf << std::endl;
    }
  } else if (algType == "test"){
    int numSample = 10;

    IncrementalInfEstimator incEstimator0;
    incEstimator0.init(&hygraph, numSample);
    std::cout << "init num inf " << incEstimator0.numInf() << std::endl;
    std::cout << "numinf for remove 82 " << incEstimator0.estimateNumInf(82) << std::endl;

    IncrementalInfEstimator incEstimator1;
    hygraph.restrictHyEdge(82);
    incEstimator1.init(&hygraph, numSample);
    std::cout << "numinf for remove 82 actual " << incEstimator1.numInf() << std::endl;

  } else {
    /* greedy algorithm */
    IncrementalInfEstimator incEstimator;
    int numSample = 10;
    incEstimator.init(&hygraph, numSample);
    for (int i = 0; i < k; ++i) {
      std::cout << "  Greedily Getting " << i << "-th Hyperedge to be restricted ..." << std::endl;
      HyEdgeID e_min; VertID minInf = VID_MAX;
//      for (HyEdgeID e = 0; e < hygraph.numHyEdges(); ++e) {
      std::vector<std::pair<HyEdgeID,uint64_t>> degrees;
      for (HyEdgeID i = 0; i < hygraph.numHyEdges(); ++i) {
        degrees.emplace_back(i, hygraph.getHyEdge(i).vertices_.size());
      }
      std::sort(degrees.begin(), degrees.end(),
                [](const std::pair<HyEdgeID,uint64_t>& left, const std::pair<HyEdgeID,uint64_t>& right){ return left.second > right.second; });
      for (int i = 0; i < 100; ++i) {
        HyEdgeID e = degrees.at(i).first;

        if (topKHyEdges.count(e) > 0) continue;
        VertID numInf = incEstimator.estimateNumInf(e);
        //std::cout << numInf << std::endl;
        if (numInf < minInf) {
          minInf = numInf;
          e_min = e;
        }
        std::cout << "\r    Finish " << e << " Hyperedges out of " << hygraph.numHyEdges() << " current min is " << e_min << ". Its inf is " << minInf << "        " << std::flush;
      }
      std::cout << std::endl;
      incEstimator.update(e_min);
      topKHyEdges.insert(e_min);
      std::cout << "  Finish. " << i << "-th Local Min is " << e_min << std::endl;
    }
  }
  std::cout << "==== Finish ====" << std::endl;
  double end = getTime();
  std::cout << "Execution Time(sec): " << (double)(end - start) << std::endl;

  /* Analyze Result */
  std::cout << "---- Start Analyze Result ----" << std::endl;
  std::cout << "!!!! Restriction Activities !!!!" << std::endl;
  std::cout << " HyEdge ID: ";
  for (auto& x: topKHyEdges) {
    std::cout << x << ",";
  }
  std::cout << std::endl;

  /* calculate # of influence (infected vertices) */
  IncrementalInfEstimator estimator;
  hygraph.restrictHyEdges(topKHyEdges);
  estimator.init(&hygraph,100);
  std::cout << "# Influence is " << estimator.numInf() << std::endl;

  std::cout << std::endl;
  std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
  std::cout << "----  Finish ----" << std::endl;

  return 0;
}