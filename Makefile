cCXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -Werror -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

define GDJDIRERR
 GDJDIR is not set at all. Please set this environment variable to point to your build - this should be either
export GDJDIR=$(PWD)
or
source setEnv.sh
if you have made appropriate changes.
For more, see README for full setup recommendations
endef

ifndef GDJDIR
$(error "$(GDJDIRERR)")	
endif

INCLUDE=-I$(GDJDIR)
LIB=-L$(GDJDIR)/lib
ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(GDJDIR)/bin
MKDIR_LIB=mkdir -p $(GDJDIR)/lib
MKDIR_OBJ=mkdir -p $(GDJDIR)/obj
MKDIR_OUTPUT=mkdir -p $(GDJDIR)/output
MKDIR_PDF=mkdir -p $(GDJDIR)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/globalDebugHandler.o obj/checkMakeDir.o obj/configParser.o obj/centralityFromInput.o lib/libATLASGDJ.so bin/gdjNTupleToHist.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirLib:
	$(MKDIR_LIB)

mkdirObj:
	$(MKDIR_OBJ)

mkdirOutput:
	$(MKDIR_OUTPUT)

mkdirPdf:
	$(MKDIR_PDF)

obj/checkMakeDir.o: src/checkMakeDir.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/checkMakeDir.C -o obj/checkMakeDir.o $(INCLUDE)

obj/globalDebugHandler.o: src/globalDebugHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/globalDebugHandler.C -o obj/globalDebugHandler.o $(ROOT) $(INCLUDE)

obj/configParser.o: src/configParser.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/configParser.C -o obj/configParser.o $(INCLUDE) $(ROOT)

obj/centralityFromInput.o: src/centralityFromInput.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o $(INCLUDE) $(ROOT)

lib/libATLASGDJ.so:
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libATLASGDJ.so obj/checkMakeDir.o obj/globalDebugHandler.o obj/configParser.o obj/centralityFromInput.o $(ROOT) $(INCLUDE)

bin/gdjNTupleToHist.exe: src/gdjNTupleToHist.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToHist.C -o bin/gdjNTupleToHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

clean:
	rm -f ./*~
	rm -f ./#*#
	rm -f bash/*~
	rm -f bash/#*#
	rm -f bin/*.exe
	rm -rf bin
	rm -f include/*~
	rm -f include/#*#
	rm -f input/*~
	rm -f input/#*#
	rm -f lib/*.so
	rm -rf lib
	rm -f obj/*.o
	rm -rf obj
	rm -f src/*~
	rm -f src/#*#
