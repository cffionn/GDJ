#!/bin/bash

lowDir=20200515
upDir=20210527
targDir=/atlasgpfs01/usatlas/data/cfmcginn/ATLASNTuples/GammaMultiJet/GDJBackup/output
mkdir -p $targDir

for i in output/*
do
    dirNum=${i#output/}
    if [[ $dirNum -lt $lowDir ]]
    then
	continue
    elif [[ $dirNum -gt $upDir ]]
    then
	continue
    fi

    cp -r $i $targDir
    rm -rf $i
#    echo $dirNum
done
