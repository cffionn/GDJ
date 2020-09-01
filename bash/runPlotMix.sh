#!/bin/bash

for i in input/mixedEvt/*.config
do
    ./bin/gdjMixedEventPlotter.exe $i
#    echo $i
done
