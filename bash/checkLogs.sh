#!/bin/bash

if [[ $# -ne 1 ]]
then
    echo "Usage: bash ./cleanLua.sh <inLogDir>. Please give input dir name as argument. exit"
    # This is janky but ensures you dont kill your terminal if you use source instead of bash
    return 1
    exit 1
else
    if [[ -d $1 ]]
    then
     dummy=0
    else
        echo "Given dir '$1' doesnt exist or is not a dir. exit"
        return 1
        exit 1
    fi
fi


for i in $1/*.log
do
    cat $i
    echo $i
    echo ""
    sleep 3
done
