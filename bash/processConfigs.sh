#!/bin/bash

dir=""
wildcard=""

if [[ $# -ge 3 ]]
then 
    echo "Greater than 2 arguments ($#) supplied."
    echo "Usage: bash processConfigs.sh <inDir> <inWildCard-optional>. return" 
    return 1 
    exit 1
elif [[ $# -eq 2 ]]
then
    dir=$1
    wildcard=$2
elif [[ $# -eq 1 ]]
then
    dir=$1
else
    echo "0 arguments supplied."
    echo "Usage: bash processConfigs.sh <inDir> <inWildCard-optional>. return" 
    return 1 
    exit 1

fi


if [[ -d $dir ]]
then
    dummy=0
else
    echo "Given '$dir' is invalid. Please provide valid directory. return"
    return 1
    exit 1
fi

searchStr=$dir/*
length=${#wildcard} 
if [[ $length -gt 0 ]]
then
    searchStr="$searchStr$wildcard"*
fi
searchStr=$searchStr.txt

while [[ $searchStr == *"//"* ]]
do
    searchStr=$(echo $searchStr | sed -e "s@//@/@g")
done

for i in $searchStr
do
    if [[ -f $i ]]
    then
	emacs -nw $i
    fi
done
