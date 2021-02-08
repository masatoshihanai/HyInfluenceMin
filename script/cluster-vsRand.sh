#!/bin/bash -x

BIN_HOME=/home/masa/workspace/cpp/HyInfluenceMin/build-cluster

GRAPH="../data/brightkite-regular.checkin ../data/gowalla-regular.checkin ../data/dataset_TSMC2014_NYC-regular.checkin ../data/dataset_TSMC2014_TKY-regular.checkin"
MEASURE="delete shrink split"
INTERVAL="weekly"

for rand in 1
do
  for interval in ${INTERVAL}
  do
    for tunit in 12 6
    do
      for graph in ${GRAPH}
      do
        for m in ${MEASURE}
        do
          for k in 10 20 30 40 50 60 70 80 90 100 1000 10000 100000
          do
            ${BIN_HOME}/HyMin -v -m ${m} -s ${rand} -r ${interval} -u ${tunit} ${graph} ${k} 2>&1 | tee -a result-new.log
          done
        done
      done
    done
  done
done
