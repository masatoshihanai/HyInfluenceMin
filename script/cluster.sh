#!/bin/bash -x

BIN_HOME=/home/masa/workspace/cpp/HyInfluenceMin/buildCluster

GRAPH="../data/brightkite.checkin ../data/gowalla.checkin ../data/dataset_TSMC2014_NYC.checkin ../data/dataset_TSMC2014_TKY.checkin"
MEASURE="delete shrink split"
INTERVAL="weekly daily"

for rand in 1
do
  for graph in ${GRAPH}
  do
    for interval in ${INTERVAL}
    do
      for tunit in 1
      do
        for m in ${MEASURE}
        do
          for k in 10 20 30 40 50 60 70 80 90 100
          do
            ${BIN_HOME}/HyMin -v -a -m ${m} -s ${rand} -r ${interval} -u ${tunit} ${graph} ${k} 2>&1 | tee -a log.log
          done
        done
      done
    done
  done
done
