#!/bin/bash

GDJDIR=/usatlas/u/cfmcginn/Projects/GDJ/
if [[ -d $GDJDIR ]]
then
    echo "GDJDIR set to '$GDJDIR'; if wrong please fix"
    export GDJDIR=$GDJDIR
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$GDJDIR/lib
else
    echo "GDJDIR given, '$GDJDIR' not found!!! Please fix" 
fi
