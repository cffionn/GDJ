#!/bin/bash

for i in `seq 0 5`
do
    ./bin/gdjToyMultiMix.exe input/toyMultiMix_10kToy$i.config
done
    
wait

for i in `seq 0 5`
do
    ./bin/gdjPlotToy.exe input/plotToy_10kToy$i.config
done
