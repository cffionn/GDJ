#!/bin/bash

inDir=/home/cfmcginn/CUBHIG/tempPlots2020/Feb21plots/

for i in $inDir/*.root
do
    cp input/ntupleToHist_TEMPLATE.txt input/ntupleToHist_TEST.txt
    sed -i -e "s@PLACEHOLDER@$i@g" input/ntupleToHist_TEST.txt
    
    echo "Running $i"
    ./bin/gdjNTupleToHist.exe input/ntupleToHist_TEST.txt
    rm input/ntupleToHist_TEST.txt
    echo ""
done
