/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#include <algorithm>
#include <deque>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <math.h>
#include <time.h>
#include "external/xoshiro256plus.hpp"

int main(int argc, char** argv) {
//  xoshiro256p::initSeed(1);
//
//  for (int i = 0; i < 100; ++i){
//    double v = xoshiro256p::to_doubleFrom0to1(xoshiro256p::next());
//    std::cout << v << std::endl;
//  }
//
//  exit(1);

  if (argc < 2) {
    std::cerr << "Data Cleaner for checkInData " << std::endl;
    std::cerr << " $ " << argv[0] << " <original checkIn-file> > <cleaned file>" << std::endl;
    std::cerr << " Example $ " << argv[0] << " ../data/loc-gowalla_totalCheckins.txt > gowalla.checkin" << std::endl;
    return 1;
  }

  std::ios::sync_with_stdio(false);


  /* Get graph data */
  char* fileName = argv[1];
  std::cout << "Open file: " << fileName << std::endl;
  std ::cout << "Read ..." << std::flush;
  std::ifstream file(fileName);
  if (file.fail()) {
    std::cerr << "!!!! File does not exist !!!! " << std::endl;
    std::cerr << "Check: " << fileName << std::endl;
    exit(1);
  }

  using VertID = long;

  using LocationID = int;
  struct Location {
    double latitude_ = 0.0;
    double longitude_ = 0.0;
  };

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
      //if (VERBOSE) std::cout << "Parsing activity category" << std::endl;
      std::getline(ss, elementStr,','); // [orig_category]
      std::getline(ss, elementStr,','); // [orig_category id]
      std::getline(ss, elementStr,','); // [category]
      std::getline(ss, elementStr,','); // [category id]
      //std::cout << elementStr << std::endl; // todo Peter
    }

    int y,mo,d,h,mi,s;
    sscanf(&checkInTimeStr.at(0), "%d-%d-%dT%d:%d:%dZ", &y, &mo, &d, &h, &mi, &s);
    tm t;
    t.tm_year = y - 1900; t.tm_mon = mo - 1; t.tm_mday = d; t.tm_hour = h; t.tm_min = mi; t.tm_sec = s; t.tm_isdst = -1;
    mktime(&t);
    checkin.checkInTime_ = t;

//    std::cout << checkInTimeStr << std::endl;
//    std::cout << "org check-in " << asctime(&t) << std::endl;
//              " year " << t.tm_year <<
//              " month " << t.tm_mon <<
//              " week " << t.tm_wday <<
//              " mday " << t.tm_mday <<
//              " yday " << t.tm_yday <<
//              " hour " << t.tm_hour << std::endl;

    checkIns.push_back(checkin);
    maxLID = maxLID < checkin.locationID_ ? checkin.locationID_: maxLID;
    maxVID = maxVID < checkin.user_ ? checkin.user_: maxVID;
  }
  file.close();
  std::cout << " ... Finish." << std::endl;
  long numCheckIn_ = checkIns.size();

  /* Output result */
  std::cout << "Analyze repeatability " << std::endl;
//  VertID prevId = 0;
//  std::unordered_map<VertID, int> map;
//  std::vector<int> freq(2100);
//  int maxFreq = 0;
//  for (auto x: checkIns) {
//    if (x.user_ == prevId) {
//      if (map.count(x.locationID_) == 0) {
//        map.emplace(x.locationID_, 1);
//      } else {
//        ++map.at(x.locationID_);
//      }
//    } else {
//      std::cout << "UserID: " << x.user_ << std::endl;
//      for (auto m: map) {
//        std::cout << "    ID " << m.first << " # " << m.second << std::endl;
//        ++freq.at(m.second);
//        if (maxFreq < m.second) maxFreq = m.second;
//      }
//      prevId = x.user_;
//      map.clear();
//    }
//  }
//  std::cout << maxFreq << std::endl;
//
//  for (auto x: freq) { std::cout << x << ","; }
//  std::cout << std::endl;

//  VertID prevId = 0;
//  std::unordered_map<LocationID, std::vector<CheckIn>> map;
//  std::vector<int> freq(1000);
//  int maxFreq = 0;
//  for (auto x: checkIns) {
//    if (x.user_ == prevId) {
//      if (map.count(x.locationID_) == 0) {
//        map.emplace(x.locationID_, std::vector<CheckIn>());
//      } else {
//        map.at(x.locationID_).push_back(x);
//      }
//    } else {
//      std::cout << "UserID: " << x.user_ << std::endl;
//      for (auto m: map) {
//        if (m.second.size() <= 1) continue;
//        std::sort(m.second.begin(), m.second.end(), [](CheckIn& r, CheckIn& l) {
//          return mktime(&r.checkInTime_) < mktime(&l.checkInTime_);});
//        int prev = mktime(&m.second.at(0).checkInTime_);
//        for (auto x: m.second) {
//          if (difftime(mktime(&x.checkInTime_),prev)/60/60/24 == 0.0) continue;
//          std::cout << "     t " << difftime(mktime(&x.checkInTime_),prev)/60/60/24 << std::endl;
//          ++freq.at((int)difftime(mktime(&x.checkInTime_),prev)/60/60/24);
//          prev = mktime(&x.checkInTime_);
//        }
//        std::cout << std::endl;
//      }
//      prevId = x.user_;
//      map.clear();
//    }
//  }
//  std::cout << maxFreq << std::endl;
//
//  for (auto x: freq) { std::cout << x << ","; }
//  std::cout << std::endl;

  VertID prevId = 0;
  std::unordered_map<LocationID, std::vector<CheckIn>> map;
  std::vector<int> freq(10000);
  int maxFreq = 0;
  int numRpEvents = 0;
  int rm = 0;
  int rm0 = 0;
  int rm1 = 0;
  int maxlength = 0;
  std::vector<int> freqWeek(7);
  std::vector<int> freqNoRep(7);
  for (auto x: checkIns) {
    if (x.user_ == prevId) {
      if (map.count(x.locationID_) == 0) {
        map.emplace(x.locationID_, std::vector<CheckIn>());
      } else {
        map.at(x.locationID_).push_back(x);
      }
    } else {
      std::cout << "UserID: " << x.user_ << std::endl;
      for (auto m: map) {
        if (m.second.size() <= 1) continue;
        ++numRpEvents;
        std::sort(m.second.begin(), m.second.end(), [](CheckIn& r, CheckIn& l) {
          return mktime(&r.checkInTime_) < mktime(&l.checkInTime_);});
        int prev = mktime(&m.second.at(0).checkInTime_);
        double prevLonti = m.second.at(0).location_.longitude_;
        double prevLati = m.second.at(0).location_.latitude_;
        int minDiff = std::numeric_limits<int>::max();
        int maxDiff = 0;
        std::vector<int> diff;
        int max = 0;
        int min = std::numeric_limits<int>::max();

        for (auto x: m.second) {
          if (difftime(mktime(&x.checkInTime_),prev)/60/60/24 == 0.0) continue;
          std::cout << " location " << x.locationID_ << " time " << asctime(&x.checkInTime_) << "         from prev " << difftime(mktime(&x.checkInTime_),prev)/60/60/24 << " dst " << sqrt((x.location_.latitude_ - prevLati)*(x.location_.latitude_ - prevLati) + (x.location_.longitude_ - prevLonti)*(x.location_.longitude_ - prevLonti)) << " coordinate: " << x.location_.latitude_ << "," << x.location_.longitude_ << std::endl;
          if ((int)difftime(mktime(&x.checkInTime_),prev)/60/60/24 > maxFreq) {
            maxFreq = (int)difftime(mktime(&x.checkInTime_),prev)/60/60/24;
          }
          if ((int)difftime(mktime(&x.checkInTime_),prev)/60/60/24 > maxDiff) {
            maxDiff = (int)difftime(mktime(&x.checkInTime_),prev)/60/60/24;
          }
          if ((int)difftime(mktime(&x.checkInTime_),prev)/60/60/24 < minDiff) {
            minDiff = (int)difftime(mktime(&x.checkInTime_),prev)/60/60/24;
          }

          if (max < mktime(&x.checkInTime_)) {
            max = mktime(&x.checkInTime_);
          }

          if (min > mktime(&x.checkInTime_)) {
            min = mktime(&x.checkInTime_);
          }

          diff.push_back((int)difftime(mktime(&x.checkInTime_),prev)/60/60/24);

          prev = mktime(&x.checkInTime_);

          ++freqWeek.at(x.checkInTime_.tm_wday);
        }
        std::cout << std::endl;
        std::cout << "max diff " << maxDiff << std::endl;
        if (diff.size() != 0) {
          std::sort(diff.begin(), diff.end());
          std::cout << "mid diff " << diff.at(diff.size()/2) << std::endl;
        }
        std::cout << "min diff " << minDiff << std::endl;

        std::cout << "max " << max << std::endl;
        std::cout << "min " << min << std::endl;

        if (minDiff > 8) ++rm;
        if (maxDiff == 0) ++rm0;
        if (diff.size() > 1 && (max-min)/60/60/24 == 0) {
          ++rm1;
          std::cout << " non rep " << std::endl;

          for (auto x: m.second) {
            ++freqNoRep.at(x.checkInTime_.tm_wday);
          }
        }
        if (maxlength < m.second.size()) {
          maxlength = m.second.size();
        }

        ++freq.at(m.second.size());

        std::cout << std::endl;
      }
      prevId = x.user_;
      map.clear();
    }
  }
  freq.resize(maxlength+1);

  std::cout << "max frequency " << maxFreq << std::endl;
  std::cout << "numRpEvents " << numRpEvents << std::endl;
  std::cout << "org events " << checkIns.size() << std::endl;
  std::cout << "rm over 7: " << rm << std::endl;
  std::cout << "rmm less than 0: " << rm0 << std::endl;
  std::cout << "rm1 less than 0 between max and min: " << rm1 << std::endl;

//  for (auto x: freq) { std::cout << x << ","; }
//  std::cout << std::endl;

  for (auto w: freqWeek) {
    std::cout << w << ",";
  }
  std::cout << std::endl;

  for (auto w: freqNoRep) {
    std::cout << w << ",";
  }
  std::cout << std::endl;
}