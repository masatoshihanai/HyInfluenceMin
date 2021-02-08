/* 
 *  Copyright (c) 2020 Masatoshi Hanai
 *
 *  This software is released under MIT License.
 *  See LICENSE.
 *
 */

#include "hygraph.hpp"
#include "external/Eigen/Core"
#include "external/Eigen/Dense"
#include "external/Eigen/Eigenvalues"
#include "wtime.hpp"

class EigenVectorCentrality {
 public:
  void getTopK(HyGraph& hygraph, HyIDSet& topK, uint64_t k) {
#ifdef HAS_OMP
    Eigen::initParallel();
#endif

    if (hygraph.numHyEdges() > 500000) {
      std::cout << "!!! Too big data for Eigen !!! Skip computation " << std::endl;
      return;
    }

    // todo fix to include repeat interval
    std::string bFileName = hygraph.fName_ + std::to_string(TIME_UNIT_HOUR) + "-eigenvec.bin";
    std::vector<std::pair<HyEdgeID, double>> idValuePairs(hygraph.numHyEdges());
    std::ifstream ifs(bFileName);
    if (ifs.is_open()) {
      std::cout << "Eigen vector exists. Read from file " << std::endl;
      ifs.read((char*)&idValuePairs[0], (sizeof(HyEdgeID) + sizeof(double) * idValuePairs.size()));
      ifs.close();
    } else {
      ifs.close();
      std::cout << "Build Matrix" << std::endl;
      Eigen::MatrixXd mmT(hygraph.numHyEdges(), hygraph.numHyEdges());
      {
        Eigen::MatrixXd m(hygraph.numHyEdges(), hygraph.numVert());
        for (uint64_t i = 0; i < hygraph.numHyEdges(); ++i) {
          for (auto v: hygraph.getHyEdge(i).vertices_) {
            m(i, v) = 1.0;
          }
        }

        auto mT = m.transpose();
        mmT = m * mT;
      }

      double start = getTime();
      uint64_t numItr = 100;

      Eigen::Vector<double, Eigen::Dynamic> v_prev = Eigen::Vector<double, Eigen::Dynamic>::Random(
        mmT.cols()).cwiseAbs();
      Eigen::Vector<double, Eigen::Dynamic> v = Eigen::Vector<double, Eigen::Dynamic>(mmT.cols());
      std::cout << "Start Power Method" << std::endl;
      for (uint64_t i = 0; i < numItr; ++i) {
        v.swap(mmT * v_prev);
        v.normalize();
        v_prev.swap(v);
        std::cout << "\r    Iteration " << i << " : " << (getTime() - start) << std::flush;
      }
      std::cout << std::endl;

      for (uint i = 0; i < idValuePairs.size(); ++i) {
        idValuePairs.at(i).first = i;
        idValuePairs.at(i).second = v(i);
      }

      std::sort(idValuePairs.begin(), idValuePairs.end(),
                [](std::pair<HyEdgeID, double> &l, std::pair<HyEdgeID, double> &r) {
                  return l.second > r.second;
                });
      double end = getTime();
      std::cout << "Execution Time for Eigen Centrality: " << end - start << std::endl;

      std::cout << "Write result to " << bFileName << std::endl;
      std::ofstream file;
      file.open(bFileName, std::ios::out | std::ios::binary);

      file.write((char*) &idValuePairs[0], (sizeof(HyEdgeID) + sizeof(double))*idValuePairs.size());

      file.close();
    }

    for (uint i = 0; i < k; ++i) {
//      std::cout << "i: " << i << " " << idValuePairs.at(i).first << " " << idValuePairs.at(i).second << std::endl;
      topK.insert(idValuePairs.at(i).first);
    }
  }
};

//int main(int argc, char** argv) {
//  int size = std::atol(argv[1]);
//  std::cout << "test EigenVector size " << size << std::endl;
//
////  Eigen::MatrixXd mmt = Eigen::MatrixXd::Random(size,size);
//
////  /* (12*7) vector */
//  Eigen::MatrixXd m = Eigen::MatrixXd(12,7);
//  m <<
//    1,1,0,0,1,0,0,
//    1,0,0,1,1,0,0,
//    0,1,1,0,1,0,0,
//    0,0,1,1,1,0,0,
//    1,1,0,0,0,1,0,
//    1,0,1,0,0,1,0,
//    1,0,0,1,0,1,0,
//    0,1,0,1,0,1,0,
//    1,0,0,1,0,0,1,
//    1,0,1,0,0,0,1,
//    0,1,0,1,0,0,1,
//    0,1,1,0,0,0,1;
//
//  auto mt = m.transpose();
//  std::cout << "m is" << std::endl;
//  std::cout << m << std::endl;
//  std::cout << std::endl;
//
//  std::cout << "transpose(m) is " << std::endl;
//  std::cout << mt << std::endl;
//  std::cout << std::endl;
//
//  std::cout << "m * transpose(m) is " << std::endl;
//  std::cout << (m*mt) << std::endl;
//  std::cout << std::endl;
//  auto mmt = m*mt;
//
//  // Exact method
//  Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> s(mmt);
//  auto eval = s.eigenvalues();
//  std::cout << "eigen value \n" << eval << std::endl;
//  std::cout << eval(0) << std::endl;
////  for (auto x: eval) {
////    std::cout << x << std::endl;
////  }
//
//  auto evec = s.eigenvectors();
//  std::cout << "eigen vector \n";// << evec << std::endl;
//
//  for (int i = 0; i < 12; ++i) {
//    std::cout << evec(i,0).real() / evec(0,0).real() << std::endl;
//  }
//
//  // Power method
//  clock_t start = clock();
//  uint64_t numItr = 100;
//
//  std::cout << "mmt size" << mmt.cols() << std::endl;
//  Eigen::Vector<double, Eigen::Dynamic> v_prev = Eigen::Vector<double, Eigen::Dynamic>::Random(mmt.cols()).cwiseAbs();
//  Eigen::Vector<double, Eigen::Dynamic> v = Eigen::Vector<double, Eigen::Dynamic>(mmt.cols());
//
//  for (uint64_t i = 0; i < numItr; ++i) {
//    v.swap(mmt*v_prev);
//    v.normalize();
//    v_prev.swap(v);
//  }
//  clock_t end = clock();
//  std::cout << "power method's eigen vector \n" << std::endl;
//  std::cout << v << std::endl;
//
//  std::vector<std::pair<HyEdgeID, double>> idValuePairs(mmt.cols());
//  for (uint i = 0; i < idValuePairs.size(); ++i) {
//    idValuePairs.at(i).first = i;
//    idValuePairs.at(i).second = v(i);
//  }
//
//  std::sort(idValuePairs.begin(), idValuePairs.end(), [](std::pair<HyEdgeID, double>& l, std::pair<HyEdgeID, double>& r) {
//    return l.second > r.second;
//  });
//
//  for (auto x: idValuePairs) {
//    std::cout << "id " << x.first << " val " << x.second << std::endl;
//  }
//
//
//  std::cout << "time to compute for " << size*size << " " << (double) (end - start)/CLOCKS_PER_SEC << std::endl;
//}
