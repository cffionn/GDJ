#!/bin/bash

#JETR,INRVAL

DATE=`date +%Y%m%d`

mkdir -p logs/$DATE/

rs=(2 4)
#files=(ntupleToHist_PbPbMC ntupleToHist_PPMC ntupleToHist_PbPbData ntupleToHist_PPData)
files=(ntupleToHist_PPMC ntupleToHist_PPData)

for i in ${rs[@]}
do
    for j in ${files[@]}
    do
	cp input/ntupleToHist/$j.config  input/ntupleToHist/"$j"_R$i.config 
	sed -i -e "s@INRVAL@$i@g" input/ntupleToHist/"$j"_R$i.config

	./bin/gdjNTupleToHist.exe input/ntupleToHist/"$j"_R$i.config &> logs/$DATE/"$j"_R"$i"_$DATE.log &
    done
    
    wait
done

echo "RUNALL.SH COMPLETE!"
