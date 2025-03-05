#!/bin/bash

includePath=$PWD # to build from other locations, just point to the GDJ directory; /home/cfm/Projects/GDJ/ in my case

mkdir -p bin
mkdir -p lib
mkdir -p obj

#### LIBS ####
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/binFlattener.C -o obj/binFlattener.o -I$includePath  `root-config --cflags --glibs`
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o -I$includePath  `root-config --cflags --glibs`
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/checkMakeDir.C -o obj/checkMakeDir.o -I$includePath
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/configParser.C -o obj/configParser.o -I$includePath `root-config --cflags --glibs`
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/globalDebugHandler.C -o obj/globalDebugHandler.o `root-config --cflags --glibs` -I$includePath 
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/keyHandler.C -o obj/keyHandler.o -I$includePath 
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/sampleHandler.C -o obj/sampleHandler.o `root-config --cflags --glibs` -I$includePath 
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -c src/mixMachine.C -o obj/mixMachine.o `root-config --cflags --glibs` -I$includePath 
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g -fPIC -shared -o lib/libATLASGDJ.so obj/binFlattener.o obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o obj/mixMachine.o `root-config --cflags --glibs` -I$includePath 


#### PLOTTER ####
g++ -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g src/gdjPlotResults.C -o bin/gdjPlotResults.exe `root-config --cflags --glibs` -I$includePath -L/home/cfm/Projects/GDJ//lib -lATLASGDJ
