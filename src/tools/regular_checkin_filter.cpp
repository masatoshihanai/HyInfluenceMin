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
#include <unordered_set>
#include <unordered_map>
#include <vector>

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
  char *fileName = argv[1];
  std::cerr << "Open file: " << fileName << std::endl;
  std::cerr << "Read ..." << std::flush;
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
    std::string str_;
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
    checkin.str_ = std::string(line);
    std::string checkInTimeStr;
    std::string elementStr;
    std::getline(ss, elementStr, ',');
    checkin.user_ = std::stol(elementStr);

    std::getline(ss, elementStr, ',');
    checkInTimeStr = elementStr;

    std::getline(ss, elementStr, ',');
    checkin.location_.latitude_ = std::stod(elementStr);

    std::getline(ss, elementStr, ',');
    checkin.location_.longitude_ = std::stod(elementStr);

    std::getline(ss, elementStr, ',');
    checkin.locationID_ = std::stol(elementStr);

    if (ss.good()) {
      //if (VERBOSE) std::cout << "Parsing activity category" << std::endl;
      std::getline(ss, elementStr, ','); // [orig_category]
      std::getline(ss, elementStr, ','); // [orig_category id]
      std::getline(ss, elementStr, ','); // [category]
      std::getline(ss, elementStr, ','); // [category id]
    }

    int y, mo, d, h, mi, s;
    sscanf(&checkInTimeStr.at(0), "%d-%d-%dT%d:%d:%dZ", &y, &mo, &d, &h, &mi, &s);
    tm t;
    t.tm_year = y - 1900;
    t.tm_mon = mo - 1;
    t.tm_mday = d;
    t.tm_hour = h;
    t.tm_min = mi;
    t.tm_sec = s;
    t.tm_isdst = -1;
    mktime(&t);
    checkin.checkInTime_ = t;

    checkIns.push_back(checkin);
    maxLID = maxLID < checkin.locationID_ ? checkin.locationID_ : maxLID;
    maxVID = maxVID < checkin.user_ ? checkin.user_ : maxVID;
  }
  file.close();
  std::cerr << " ... Finish." << std::endl;
  long numCheckIn_ = checkIns.size();

  /* Output result */
  std::cerr << "Filter non repeat and repeat within a day." << std::endl;

  long numUser = 0;
  VertID prevId = 0;
  long regular = 0;
  long single = 0;
  long nonRep = 0;
  std::unordered_map<LocationID, std::vector<CheckIn>> gatherByLocation;
  std::vector<long> weeklyFreqNonRep(7,0);
  std::vector<long> weeklyFreqRegular(7,0);
  std::vector<long> weeklyFreqSingle(7,0);
  for (auto x: checkIns) {
    if (x.user_ == prevId) {
      // put to the same location
      if (gatherByLocation.count(x.locationID_) == 0) {
        gatherByLocation.emplace(x.locationID_, std::vector<CheckIn>());
      }
      gatherByLocation.at(x.locationID_).push_back(x);
    } else {
      ++numUser;
      prevId = x.user_;
      for (auto chins: gatherByLocation) {
        if (chins.second.size() <= 1) {
          ++weeklyFreqNonRep.at(chins.second.front().checkInTime_.tm_wday);
          ++nonRep;
          continue; // skip non repeated check-in
        }

        // sort by time
        std::sort(chins.second.begin(), chins.second.end(), [](CheckIn& r, CheckIn& l) {
          return mktime(&r.checkInTime_) < mktime(&l.checkInTime_);});

        if (difftime(mktime(&chins.second.back().checkInTime_), mktime(&chins.second.front().checkInTime_))/60/60/24 < 1) {
          for (auto x: chins.second) {
            ++weeklyFreqSingle.at(x.checkInTime_.tm_wday);
            ++single;
          }
          continue;
        };

        for (auto result: chins.second) {
          ++regular;
          ++weeklyFreqRegular.at(result.checkInTime_.tm_wday);
          std::cout << result.str_ << std::endl;
        }
      }
      gatherByLocation.clear();

      gatherByLocation.emplace(x.locationID_, std::vector<CheckIn>());
      gatherByLocation.at(x.locationID_).push_back(x);
    }
  }

  std::cerr << "regular " << regular << " single " << single << " non rep " << nonRep << " origin " << checkIns.size() << std::endl;

  std::cerr << "statistical result" << std::endl;
  std::cerr << "Reqular,";
  for (auto r: weeklyFreqRegular) { std::cerr << r << ","; }
  std::cerr << std::endl;

  std::cerr << "NonRep,";
  for (auto r: weeklyFreqNonRep) { std::cerr << r << ","; }
  std::cerr << std::endl;

  std::cerr << "Single,";
  for (auto r: weeklyFreqSingle) { std::cerr << r << ","; }
  std::cerr << std::endl;
}