#!/bin/bash

GDJDIR=/home/cfm/Projects/GDJ/
ROOUNFOLDDIR=/home/cfm/Packages/RooUnfold/RooUnfold-build/
#HEPMCLIB=/home/cfm/Packages/HepMC2/hepmc-install/lib/

if [[ $HOME == "/usatlas/u/cfmcginn" ]]
then
    HEPMCPATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_97a/HepMC/2.06.11/x86_64-centos7-gcc8-opt/lib
    if [[ $LD_LIBRARY_PATH == *"$HEPMCPATH"* ]]
    then
	dummy=0
    else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/sft.cern.ch/lcg/releases/LCG_97a/HepMC/2.06.11/x86_64-centos7-gcc8-opt/lib
    fi

    GDJDIR=/usatlas/u/cfmcginn/Projects/GDJ/
    ROOUNFOLDDIR=/usatlas/u/cfmcginn/Packages/RooUnfold-build/
fi

GDJDIRLIB="$GDJDIR"lib

export DOGLOBALDEBUGROOT=0 #from command line, initiating


if [[ -d $GDJDIR ]]
then
    echo "GDJDIR set to '$GDJDIR'; if wrong please fix"
    export GDJDIR=$GDJDIR

    
    if [[ $LD_LIBRARY_PATH == *"$GDJDIRLIB"* ]]
    then
	dummy=0
    else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GDJDIRLIB
    fi

    if [[ -d $ROOUNFOLDDIR ]]
    then
	export ROOUNFOLDDIR=$ROOUNFOLDDIR
	
	if [[ $LD_LIBRARY_PATH == *"$ROOUNFOLDDIR"* ]]
	then
	    dummy=0
	else
	    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOUNFOLDDIR
	fi
    else
	echo "NOTE: ROOUNFOLDDIR MISSING"	
    fi    
else
    echo "GDJDIR given, '$GDJDIR' not found!!! Please fix" 
fi
