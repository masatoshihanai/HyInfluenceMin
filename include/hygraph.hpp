/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#ifndef HYMINSOLVER_HYGRAPH_HPP
#define HYMINSOLVER_HYGRAPH_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <unordered_set>

//int TIME_UNIT = 86400; // 1 day = 60 * 60 * 24 = 86400 sec
//int TIME_UNIT = 86400*7*4; // 7 day = 60 * 60 * 24 * 7= 86400 sec
int TIME_UNIT = std::numeric_limits<int>::max(); // No time unit

using VertID = int;
VertID VID_MAX = std::numeric_limits<VertID>::max();
using Vertices = std::vector<VertID>;

using LocationID = int;
struct Location {
  double latitude_ = 0.0;
  double longitude_ = 0.0;
};
using Locations = std::vector<Location>;

using HyEdgeID = int;
struct HyEdge{
  double IF_;
  LocationID lID_ = 0;
  std::unordered_set<VertID> vertices_;
};
using HyEdges = std::vector<HyEdge>;

class HyGraph {
  VertID numVert_;
  std::vector<double> vertIFs_;
  std::vector<std::vector<HyEdgeID>> neighbors_;
  HyEdges hyEdges_;
  Locations locations_;
  uint64_t numCheckIn_;
  double maxMultiInfFactor;

 public:
  std::vector<HyEdgeID>& neighbors(VertID id) {
    return neighbors_.at(id);
  }

  HyEdge& getHyEdge(HyEdgeID id) {
    return hyEdges_.at(id);
  }

  double getVertIF(VertID v) {
    return vertIFs_.at(v);
  }

  double getIF(VertID v, HyEdgeID e) {
    return (double) getVertIF(v)*getHyEdge(e).IF_ / maxMultiInfFactor;
  }

  uint64_t numVert() { return numVert_; }
  uint64_t numHyEdges() { return hyEdges_.size(); }
  uint64_t numCheckIn() { return numCheckIn_; }

  void restrictHyEdge(HyEdgeID hyedge) {
    hyEdges_.at(hyedge).vertices_.clear();
  }

  void restrictHyEdges(std::unordered_set<HyEdgeID>& hyedges) {
    for (auto hyedgeID: hyedges) {
      restrictHyEdge(hyedgeID);
    }
  }

  void init(char* fileName) {
    std::ios::sync_with_stdio(false);
    /* Get graph data */
    std::cout << "Open file: " << fileName << std::endl;
    std::cout << "Read ..." << std::flush;
    std::ifstream file(fileName);
    if (file.fail()) {
      std::cerr << "!!!! File does not exist !!!! " << std::endl;
      std::cerr << "Check: " << fileName << std::endl;
      exit(1);
    }

    struct CheckIn {
      VertID user_;
      LocationID locationID_;
      Location location_;
      time_t checkInTime_;
    };
    std::deque<CheckIn> checkIns;

    std::string line;
    uint64_t count = -1;
    LocationID maxLID = 0;
    VertID maxVID = 0;
    while (std::getline(file, line)) {
      ++count;
      if (line[0] == '#') continue;
      if (line[0] == '%') continue;
      if (line[0] == '/') continue;

      std::istringstream iss(line);
      CheckIn checkin;
      std::string checkInTimeStr;
      if (!(iss >> checkin.user_
                >> checkInTimeStr
                >> checkin.location_.latitude_
                >> checkin.location_.longitude_
                >> checkin.locationID_)) {
        std::cerr << "Fail to read files. Check line " << count << std::endl;
        std::cerr << "Clean up File first. " << std::endl;
        exit(1);
      }
      int y,mo,d,h,mi,s;
      sscanf(&checkInTimeStr.at(0), "%d-%d-%dT%d:%d:%dZ", &y, &mo, &d, &h, &mi, &s);
      tm t;
      t.tm_year = y; t.tm_mon = mo; t.tm_mday = d; t.tm_mday = d; t.tm_hour = h; t.tm_min = mi; t.tm_sec = s; t.tm_isdst = -1;
      checkin.checkInTime_ = mktime(&t);
      checkIns.push_back(checkin);
      maxLID = maxLID < checkin.locationID_ ? checkin.locationID_: maxLID;
      maxVID = maxVID < checkin.user_ ? checkin.user_: maxVID;
    }
    file.close();
    std::cout << " ... Finish." << std::endl;
    numCheckIn_ = checkIns.size();

    /* Construct Graph */
    numVert_ = maxVID + 1;
    locations_.resize(maxLID+1);
    for (auto x: checkIns) {
      locations_.at(x.locationID_) = x.location_;
    }

    hyEdges_.reserve(checkIns.size());
    std::vector<std::unordered_set<HyEdgeID>> location2HyEdge(maxLID+1);
    std::deque<time_t> hyEdge2time;
    int i = 0;
    for (auto x: checkIns) {
      if (location2HyEdge.at(x.locationID_).size() == 0) {
        location2HyEdge.at(x.locationID_).insert(i);
        HyEdge hyEdge;
        hyEdge.lID_ = x.locationID_;
        hyEdge.vertices_.insert(x.user_);
        hyEdges_.push_back(hyEdge);
        hyEdge2time.push_back(x.checkInTime_ - (x.checkInTime_ % TIME_UNIT));
        ++i;
      } else {
        bool assigned = false;
        for (auto hyedgeid: location2HyEdge.at(x.locationID_)) {
          if (x.checkInTime_ > hyEdge2time.at(hyedgeid) && x.checkInTime_ - hyEdge2time.at(hyedgeid) < TIME_UNIT) {
            HyEdge &hyEdge = hyEdges_.at(hyedgeid);
            hyEdge.vertices_.insert(x.user_);
            assigned = true;
            break;
          }
        }
        if (!assigned) {
          location2HyEdge.at(x.locationID_).insert(i);
          HyEdge hyEdge;
          hyEdge.lID_ = x.locationID_;
          hyEdge.vertices_.insert(x.user_);
          hyEdges_.push_back(hyEdge);
          hyEdge2time.push_back(x.checkInTime_ - (x.checkInTime_ % TIME_UNIT));
          ++i;
        }
      }
    }

    neighbors_.resize(numVert_);
    for (HyEdgeID e = 0; e < hyEdges_.size(); ++e) {
      for (auto& v: hyEdges_.at(e).vertices_) {
        neighbors_.at(v).push_back(e);
      }
    }

    /* init infection factor */
    // todo init infection factor from file
    vertIFs_ = std::vector<double>(numVert_);
    std::mt19937 randgen;
    randgen.seed(99);
    std::lognormal_distribution<float> lognorm(0.0, 0.25);
    //xoshiro256p::initSeed(99);
    for(auto& x: vertIFs_) {
      //x = 1.0;//lognorm(randgen);//
      x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
    }
    for (auto& x: hyEdges_) {
      //x.IF_ = 1.0;//lognorm(randgen);
      x.IF_ = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
    }

    /* normalize infection factor by the max-degree vertex */
    maxMultiInfFactor = 0;
    for (VertID v = 0; v < numVert(); ++v) {
      double infFactor = 0;
      for (HyEdgeID i = 0; i < neighbors(v).size(); ++i) {
        infFactor += hyEdges_.at(neighbors(v).at(i)).IF_*vertIFs_.at(v);
      }
      if (infFactor > maxMultiInfFactor) {
        maxMultiInfFactor = infFactor;
      }
    }
    if (maxMultiInfFactor < 1) maxMultiInfFactor = 1;
    std::cout << "maxMultiInfFactor" << maxMultiInfFactor << std::endl;
    //maxMultiInfFactor = maxMultiInfFactor;

    // todo remove
    int max = 0;
    for (auto& x: hyEdges_) {
      if (x.vertices_.size() > max) {
        max = x.vertices_.size();
      }
    }
    std::vector<int> sizeDistribution(max+1, 0);
    for (auto& x: hyEdges_) {
      ++sizeDistribution.at(x.vertices_.size());
    }
    for (auto& x: sizeDistribution) {
      std::cout << x << ",";
    }
    std::cout << std::endl;

    max = 0;
    std::vector<int> degree(numVert(), 0);
    for (VertID v = 0; v < numVert(); ++v) {
      for (HyEdgeID i = 0; i < neighbors(v).size(); ++i) {
        ++degree.at(v);
      }
      if (degree.at(v) > max) {
        max = degree.at(v);
      }
    }

    std::vector<int> degreeDistribution(max+1, 0);
    for (VertID v = 0; v < numVert(); ++v) {
      ++degreeDistribution.at(degree.at(v));
    }
    for (auto& x: degreeDistribution) {
      std::cout << x << ",";
    }
    std::cout << std::endl;
    std::vector<int> actDegreeDistribution(max+1, 0);


  }
};


#endif //HYMINSOLVER_HYGRAPH_HPP
