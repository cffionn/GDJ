#!/bin/bash

files=(plotRes_R2_DATA_TEST.config plotRes_R4_DATA_TEST.config plotRes_R2_MC_TEST.config plotRes_R4_MC_TEST.config)

for i in "${files[@]}"
do
    ./bin/gdjResponsePlotter.exe input/plotRes/$i
done
	 
