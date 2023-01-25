#!/bin/bash

GDJDIR=/home/cfmcginn/Projects/GDJ/
GDJDIRLIB="$GDJDIR"lib
ROOUNFOLDDIR=/home/cfmcginn/Packages/RooUnfold/RooUnfold-build/
HEPMCLIB=/home/cfmcginn/Packages/HepMC2/hepmc-install/lib/
#THIS IS A DUMMY FOR MAKE - UNFOLD HAPPENS OFFLINE
#ROOUNFOLDDIR=$GDJDIR

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

    if [[ $LD_LIBRARY_PATH == *"$HEPMCLIB"* ]]
    then
	dummy=0
    else
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMCLIB
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

#For debugging
export DOGLOBALDEBUGROOT=0
