/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

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
  std::ifstream file(fileName);
  if (file.fail()) {
    std::cerr << "!!!! File does not exist !!!! " << std::endl;
    std::cerr << "Check " << fileName << std::endl;
    exit(1);
  }

  struct Checkin {
    std::string user;
    std::string checkInTime;
    double latitude, lontitude;
    std::string locationID;
  };
  std::deque<Checkin> checkins;
  std::unordered_map<std::string, uint64_t> userSet;
  std::unordered_map<std::string, uint64_t> locationSet;
  uint64_t count = -1;
  std::string line;
  uint64_t numSkip = 0;
  while (std::getline(file, line)) {
    ++count;
    if (line[0] == '#') continue;
    if (line[0] == '%') continue;
    if (line[0] == '/') continue;

    std::istringstream iss(line);
    Checkin checkin;
    std::string locationID;
    if (!(iss
      >> checkin.user
      >> checkin.checkInTime
      >> checkin.latitude
      >> checkin.lontitude
      >> checkin.locationID)) {
      std::cerr << "Skip to read line " << count << " \"" << line << "\""<< std::endl;
      ++numSkip;
      continue;
    } else if (checkin.latitude == 0 || checkin.lontitude == 0) {
      std::cerr << "Skip to read line " << count << " \"" << line << "\""<< std::endl;
      ++numSkip;
      continue;
    }

    if (checkin.checkInTime.size() != 20) {
      std::cout << "Error to parse time and date at line " << count << " \"" << line << "\""<< std::endl;
      exit(1);
    }

    checkins.push_back(checkin);
    if (userSet.count(checkin.user) == 0) {
      userSet.emplace(checkin.user, userSet.size());
    }
    if (locationSet.count(checkin.locationID) == 0) {
      locationSet.emplace(checkin.locationID, locationSet.size());
    }
  }
  file.close();

  /* Output result */
  std::cerr << "# CheckIn: " << checkins.size() << " (Before Clean up) " << count << ". # Users: "<< userSet.size() << ". # Locations: " << locationSet.size() << std::endl;
  std::cout << "## # CheckIn: " << checkins.size() << " (Before Clean up) " << count << ". # Users: "<< userSet.size() << ". # Locations: " << locationSet.size() << std::endl;
  std::cout << "## [user] [check-in time(20YY-MM-DD-HH:MM:SS)] [latitude] [longitude] [location id]" << std::endl;
  for (auto x: checkins) {
    std::cout << userSet.at(x.user) << " " << x.checkInTime  << " " << x.latitude << " " << x.lontitude << " " << locationSet.at(x.locationID) << std::endl;
  }
  std::cerr << "Finish" << std::endl;
}