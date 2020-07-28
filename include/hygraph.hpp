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
#include <string>
#include <time.h>
#include <unordered_set>
#include <unordered_map>

int TIME_UNIT_HOUR = 12; // Chech-ins within 12 hours are considered to be the same activity

using VertID = int;
VertID VID_MAX = std::numeric_limits<VertID>::max();
using Vertices = std::vector<VertID>;
using VertSet = std::unordered_set<VertID>;

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
  //std::unordered_set<VertID> vertices_;
  Vertices vertices_;
};
using HyEdges = std::vector<HyEdge>;
using HyIDSet = std::unordered_set<HyEdgeID>;

class HyGraph {
  VertID numVert_;
  std::vector<double> vertIFs_;
  std::vector<std::vector<HyEdgeID>> neighbors_;
  HyEdges hyEdges_;
  Locations locations_;
  uint64_t numCheckIn_;
  std::vector<double> maxMultiInfFactor;

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

  double getIF(VertID v, VertID u, HyEdgeID e) {
    return (double) getVertIF(v)*getHyEdge(e).IF_ / maxMultiInfFactor.at(u);
  }

  double genRandThr() {
    return xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
  }

 private:
  std::vector<std::unordered_map<VertID, int>> restGroup_;
 public:
  bool isRestricted(HyEdgeID edge, VertID src, VertID dst) {
    return restGroup_.at(edge).at(src) == restGroup_.at(edge).at(dst);
  }

  uint64_t numVert() { return numVert_; }
  uint64_t numHyEdges() { return hyEdges_.size(); }
  uint64_t numCheckIn() { return numCheckIn_; }

  void init(char* fileName, const std::string& restType, const std::string& repeatInterval) {
    std::ios::sync_with_stdio(false);
    /* Get graph data */
    std::cout << "Open file: " << fileName << std::endl;
    std ::cout << "Read ..." << std::flush;
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
      tm checkInTime_;
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

      std::stringstream ss(line);
      CheckIn checkin;
      std::string checkInTimeStr;
      std::string elementStr;
      std::getline(ss, elementStr,',');
      checkin.user_ = std::stol(elementStr);

      std::getline(ss, elementStr,',');
      checkInTimeStr = elementStr;

      std::getline(ss, elementStr,',');
      checkin.location_.latitude_ = std::stod(elementStr);

      std::getline(ss, elementStr,',');
      checkin.location_.longitude_ = std::stod(elementStr);

      std::getline(ss, elementStr,',');
      checkin.locationID_ = std::stol(elementStr);

      if (ss.good()) {
        if (VERBOSE) std::cout << "Parsing activity category" << std::endl;
        std::getline(ss, elementStr,','); // [orig_category]
        std::getline(ss, elementStr,','); // [orig_category id]
        std::getline(ss, elementStr,','); // [category]
        std::getline(ss, elementStr,','); // [category id]
        std::cout << elementStr << std::endl; // todo Peter
      }

      int y,mo,d,h,mi,s;
      sscanf(&checkInTimeStr.at(0), "%d-%d-%dT%d:%d:%dZ", &y, &mo, &d, &h, &mi, &s);
      tm t;
      t.tm_year = y; t.tm_mon = mo; t.tm_mday = d; t.tm_hour = h; t.tm_min = mi; t.tm_sec = s; t.tm_isdst = -1;
      mktime(&t);
      checkin.checkInTime_ = t;
//      std::cout << "org check-in " <<
//                " year " << t.tm_year <<
//                " month " << t.tm_mon <<
//                " week " << t.tm_wday <<
//                " mday " << t.tm_mday <<
//                " yday " << t.tm_yday <<
//                " hour " << t.tm_hour << std::endl;
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
    std::vector<std::unordered_set<VertID>> hyEdge2Verts;
    hyEdge2Verts.reserve(checkIns.size());

    std::vector<std::unordered_set<HyEdgeID>> location2HyEdge(maxLID+1);
    std::deque<tm> hyEdge2time;
    int i = 0;
    bool isWeekly = repeatInterval == "weekly" ? true: false;
    for (auto x: checkIns) {
      if (location2HyEdge.at(x.locationID_).size() == 0) {
        location2HyEdge.at(x.locationID_).insert(i);
        HyEdge hyEdge;
        hyEdge.lID_ = x.locationID_;
        hyEdges_.push_back(hyEdge);
        std::unordered_set<HyEdgeID> sets;
        sets.insert(x.user_);
        hyEdge2Verts.push_back(sets);
        hyEdge2time.push_back(x.checkInTime_);
        ++i;
      } else {
        bool assigned = false;
        for (auto hyedgeid: location2HyEdge.at(x.locationID_)) {
          if (isWeekly) {
            /* weekly activity */
//            std::cout << "weekly check-in " <<
//              " year " << x.checkInTime_.tm_year <<
//              " month " << x.checkInTime_.tm_mon <<
//              " week " << x.checkInTime_.tm_wday <<
//              " mday " << x.checkInTime_.tm_mday <<
//              " yday " << x.checkInTime_.tm_yday <<
//              " hour " << x.checkInTime_.tm_hour << std::endl;
            int week = hyEdge2time.at(hyedgeid).tm_wday;
            int hour = hyEdge2time.at(hyedgeid).tm_hour / TIME_UNIT_HOUR;
            int newWeek = x.checkInTime_.tm_wday;
            int newHour = x.checkInTime_.tm_hour / TIME_UNIT_HOUR;
            if (week == newWeek && hour == newHour) {
              hyEdge2Verts.at(hyedgeid).insert(x.user_);
              assigned = true;
              break;
            }
          } else {
//            std::cout << "daily check-in " <<
//                      " year " << x.checkInTime_.tm_year <<
//                      " month " << x.checkInTime_.tm_mon <<
//                      " week " << x.checkInTime_.tm_wday <<
//                      " mday " << x.checkInTime_.tm_mday <<
//                      " yday " << x.checkInTime_.tm_yday <<
//                      " hour " << x.checkInTime_.tm_hour << std::endl;
            /* daily activity */
            int hour = hyEdge2time.at(hyedgeid).tm_hour / TIME_UNIT_HOUR;
            int newHour = x.checkInTime_.tm_hour / TIME_UNIT_HOUR;
            if (hour == newHour) {
              hyEdge2Verts.at(hyedgeid).insert(x.user_);
              assigned = true;
              break;
            }
          }
        }
        if (!assigned) {
          location2HyEdge.at(x.locationID_).insert(i);
          HyEdge hyEdge;
          hyEdge.lID_ = x.locationID_;
          hyEdges_.push_back(hyEdge);
          std::unordered_set<HyEdgeID> sets;
          sets.insert(x.user_);
          hyEdge2Verts.push_back(sets);
          hyEdge2time.push_back(x.checkInTime_);
          ++i;
        }
      }
    }

    for (uint64_t i = 0; i < hyEdge2Verts.size(); ++i) {
      for (auto v: hyEdge2Verts.at(i)) {
        hyEdges_.at(i).vertices_.push_back(v);
      }
    }
    hyEdge2Verts.clear();

    // todo
    int restrictionType;
    const int DELETE = 1; const int SHRINK = 2; const int SPLIT = 3;
    if (restType == "delete") {
      restrictionType = DELETE;
    } else if (restType == "shrink") {
      restrictionType = SHRINK;
    } else if (restType == "split") {
      restrictionType = SPLIT;
    } else {
      std::cout << restType << " is invalid measure type. set delete " << std::endl;
      restrictionType = DELETE;
    }
    restGroup_.resize(hyEdges_.size());
    for (uint64_t i = 0; i < restGroup_.size(); ++i) {
      for (uint64_t j = 0; j < hyEdges_.at(i).vertices_.size(); ++j) {
        VertID v = hyEdges_.at(i).vertices_.at(j);
        if (restrictionType == DELETE) {
          restGroup_.at(i).emplace(v, 0);
        } else if (restrictionType == SHRINK) {
          // todo shrink to half size
          if (j % 2) {
            restGroup_.at(i).emplace(v, 0);
          } else {
            restGroup_.at(i).emplace(v, j);
          }
        } else if (restrictionType == SPLIT) {
          // todo split to two activities
          if (j % 2) {
            restGroup_.at(i).emplace(v, 0);
          } else {
            restGroup_.at(i).emplace(v, 1);
          }
        } else {
          /* never reached */
          std::cerr << "Error at Line " << __LINE__ << " File " << __FILE__ << std::endl;
          exit(1);
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
    double min = 1.0;
    double max = 1000;
    double alpha = 2.1;
    //xoshiro256p::initSeed(99);
    auto zipf = [&]{ return std::pow((std::pow(max, (1.0-alpha)) - std::pow(min, (1.0-alpha))) * xoshiro256p::to_doubleFrom0to1(xoshiro256p::next()) + std::pow(min,(1.0-alpha)), 1.0/(1.0-alpha));};

    vertIFs_ = std::vector<double>(numVert_);
    for(auto& x: vertIFs_) {
      x = 1.0;
      //x = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
      //x = zipf();
    }
    for (auto& x: hyEdges_) {
      //x.IF_ = 1;
      x.IF_ = zipf();
      //x.IF_ = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
    }

    /* normalize infection factor */
    double maxIF = 0.0;
    uint64_t maxDegree = 0;
    maxMultiInfFactor.resize(numVert());
    for (VertID v = 0; v < numVert(); ++v) {
      double infFactor = 0;
      uint64_t degree = 0;
      for (auto hyID: neighbors(v)) {
        for (auto u :hyEdges_.at(hyID).vertices_) {
          infFactor += hyEdges_.at(hyID).IF_*vertIFs_.at(u);
          ++degree;
        }
      }
      maxMultiInfFactor.at(v) = infFactor;
      if (infFactor > maxIF) {
        maxIF = infFactor;
      }
      if (degree > maxDegree) {
        maxDegree = degree;
      }
    }
//    std::cout << "maxMultiInfFactor" << (uint64_t) maxIF << " max degree " << maxDegree << std::endl;

    // todo
    if (false) {
      std::cout << "Normalized by max infection factor " << std::endl;
      for (VertID v = 0; v < numVert(); ++v) {
        maxMultiInfFactor.at(v) = maxIF;
      }
    }

    // todo remove
//    int max = 0;
//    for (auto& x: hyEdges_) {
//      if (x.vertices_.size() > max) {
//        max = x.vertices_.size();
//      }
//    }
//    std::vector<int> sizeDistribution(max+1, 0);
//    for (auto& x: hyEdges_) {
//      ++sizeDistribution.at(x.vertices_.size());
//    }
//    for (auto& x: sizeDistribution) {
//      std::cout << x << ",";
//    }
//    std::cout << std::endl;
//
//    max = 0;
//    std::vector<int> degree(numVert(), 0);
//    for (VertID v = 0; v < numVert(); ++v) {
//      for (HyEdgeID i = 0; i < neighbors(v).size(); ++i) {
//        ++degree.at(v);
//      }
//      if (degree.at(v) > max) {
//        max = degree.at(v);
//      }
//    }
//
//    std::vector<int> degreeDistribution(max+1, 0);
//    for (VertID v = 0; v < numVert(); ++v) {
//      ++degreeDistribution.at(degree.at(v));
//    }
//    for (auto& x: degreeDistribution) {
//      std::cout << x << ",";
//    }
//    std::cout << std::endl;
//    std::vector<int> actDegreeDistribution(max+1, 0);
  }
};


#endif //HYMINSOLVER_HYGRAPH_HPP
