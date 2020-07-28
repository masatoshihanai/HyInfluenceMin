#!/bin/bash -x

BIN_HOME=/home/masa/workspace/cpp/HyInfluenceMin/buildCluster

GRAPH="../data/brightkite.checkin ../data/gowalla.checkin"
MEASURE="delete shrink split"
INTERVAL="weekly daily"

for rand in 1 2 3 4 5
do
  for interval in ${INTERVAL}
  do
    for tunit in 24 12 6
    do
      for graph in ${GRAPH}
      do
        for m in ${MEASURE}
        do
          for k in 10 100 1000 10000 100000
          do
            ${BIN_HOME}/HyMin -v -m ${m} -s ${rand} -r ${interval} -u ${tunit} ${graph} ${k} 2>&1 | tee -a log.log
          done
        done
      done
    done
  done
done