#!/bin/bash

DATE=`date +%Y%m%d`
mkdir -p logs/$DATE/

rs=(2 4)
files=(dataMCRaw dataMCRaw_PP)
inDate=20200730

for i in ${rs[@]}
do
    for j in ${files[@]}
    do
	sysStr="PbPb"
	if [[ $j == *"PP"* ]]
	then
	    sysStr="PP"
	fi

	cp input/dataMCRaw/$j.config input/dataMCRaw/"$j"_R$i.config
	dataFile=output/$inDate/ntuplePreProc_"$sysStr"Data_R"$i"_HIST_$inDate.root
	mcFile=output/$inDate/ntupleMCPreProc_"$sysStr"MC_R"$i"_HIST_$inDate.root
	sed -i -e "s@INPUTDATAFILE@$dataFile@g" input/dataMCRaw/"$j"_R$i.config
	sed -i -e "s@INPUTMCFILE@$mcFile@g" input/dataMCRaw/"$j"_R$i.config
	
	./bin/gdjDataMCRawPlotter.exe input/dataMCRaw/"$j"_R$i.config &> logs/$DATE/"$j"_R"$i"_$DATE.log &
    done
done

wait

echo "RUNPLOTALL.SH COMPLETE!"
