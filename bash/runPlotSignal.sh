#!/bin/bash

DATE=`date +%Y%m%d`

mkdir -p logs/$DATE
#echo $DATE

./bin/gdjPlotSignalHist.exe input/plotSignalHist_R2.config &> logs/$DATE/signalPlot_R2_$DATE.log &
./bin/gdjPlotSignalHist.exe input/plotSignalHist_R4.config &> logs/$DATE/signalPlot_R2_$DATE.log &

wait

echo "runPlotSignal.sh COMPLETE"
