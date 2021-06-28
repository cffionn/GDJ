CXX = g++
#O3 for max optimization (go to 0 for debug)
CXXFLAGS = -Wall -O3 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11 -g
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

define ROOUNFOLDDIRERR
 ROOUNFOLDDIR is not set at all. Please set this environment variable to point to your RooUnfold - source setEnv.sh if you have made appropriate changes. For more, see README for full setup recommendations
endef

ifndef GDJDIR
$(error "$(GDJDIRERR)")	
endif

ifndef ROOUNFOLDDIR
$(error "$(ROOUNFOLDDIRERR)")	
endif

ROOUNFOLDDIR=/home/cfmcginn/Packages/RooUnfold/RooUnfold-build/
INCLUDE=-I$(GDJDIR) -I$(ROOUNFOLDDIR)
LIB=-L$(GDJDIR)/lib
ROOUNFOLDLIB=-L$(ROOUNFOLDDIR) -lRooUnfold

ROOT=`root-config --cflags --glibs`

MKDIR_BIN=mkdir -p $(GDJDIR)/bin
MKDIR_LIB=mkdir -p $(GDJDIR)/lib
MKDIR_OBJ=mkdir -p $(GDJDIR)/obj
MKDIR_OUTPUT=mkdir -p $(GDJDIR)/output
MKDIR_PDF=mkdir -p $(GDJDIR)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o obj/mixMachine.o lib/libATLASGDJ.so bin/gdjNtuplePreProc.exe bin/gdjNTupleToHist.exe bin/gdjNTupleToMBHist.exe bin/gdjHistDumper.exe bin/gdjGammaJetResponsePlot.exe bin/gdjMixedEventPlotter.exe bin/gdjControlPlotter.exe bin/gdjResponsePlotter.exe bin/gdjDataMCRawPlotter.exe bin/grlToTex.exe bin/testKeyHandler.exe bin/testSampleHandler.exe bin/gdjPlotMBHist.exe bin/gdjHistToUnfold.exe bin/gdjPlotUnfoldDiagnostics.exe bin/gdjPlotResults.exe
#bin/gdjNTupleToSignalHist.exe bin/gdjPlotSignalHist.exe bin/gdjToyMultiMix.exe bin/gdjPlotToy.exe

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

obj/centralityFromInput.o: src/centralityFromInput.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/centralityFromInput.C -o obj/centralityFromInput.o $(INCLUDE) $(ROOT)

obj/checkMakeDir.o: src/checkMakeDir.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/checkMakeDir.C -o obj/checkMakeDir.o $(INCLUDE)

obj/configParser.o: src/configParser.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/configParser.C -o obj/configParser.o $(INCLUDE) $(ROOT)

obj/globalDebugHandler.o: src/globalDebugHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/globalDebugHandler.C -o obj/globalDebugHandler.o $(ROOT) $(INCLUDE)

obj/keyHandler.o: src/keyHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/keyHandler.C -o obj/keyHandler.o $(INCLUDE)

obj/sampleHandler.o: src/sampleHandler.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/sampleHandler.C -o obj/sampleHandler.o $(ROOT) $(INCLUDE)

obj/mixMachine.o: src/mixMachine.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/mixMachine.C -o obj/mixMachine.o $(ROOT) $(INCLUDE)

lib/libATLASGDJ.so:
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libATLASGDJ.so obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o obj/mixMachine.o $(ROOT) $(INCLUDE)

bin/gdjNtuplePreProc.exe: src/gdjNtuplePreProc.C
	$(CXX) $(CXXFLAGS) src/gdjNtuplePreProc.C -o bin/gdjNtuplePreProc.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNTupleToHist.exe: src/gdjNTupleToHist.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToHist.C -o bin/gdjNTupleToHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNTupleToSignalHist.exe: src/gdjNTupleToSignalHist.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToSignalHist.C -o bin/gdjNTupleToSignalHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjNTupleToMBHist.exe: src/gdjNTupleToMBHist.C
	$(CXX) $(CXXFLAGS) src/gdjNTupleToMBHist.C -o bin/gdjNTupleToMBHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

#bin/quickEventIso.exe: src/quickEventIso.C
#	$(CXX) $(CXXFLAGS) src/quickEventIso.C -o bin/quickEventIso.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

#bin/getEntryInAOD.exe: src/getEntryInAOD.C
#	$(CXX) $(CXXFLAGS) src/getEntryInAOD.C -o bin/getEntryInAOD.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotMBHist.exe: src/gdjPlotMBHist.C
	$(CXX) $(CXXFLAGS) src/gdjPlotMBHist.C -o bin/gdjPlotMBHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotSignalHist.exe: src/gdjPlotSignalHist.C
	$(CXX) $(CXXFLAGS) src/gdjPlotSignalHist.C -o bin/gdjPlotSignalHist.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHistDumper.exe: src/gdjHistDumper.C
	$(CXX) $(CXXFLAGS) src/gdjHistDumper.C -o bin/gdjHistDumper.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjGammaJetResponsePlot.exe: src/gdjGammaJetResponsePlot.C
	$(CXX) $(CXXFLAGS) src/gdjGammaJetResponsePlot.C -o bin/gdjGammaJetResponsePlot.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjMixedEventPlotter.exe: src/gdjMixedEventPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjMixedEventPlotter.C -o bin/gdjMixedEventPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjControlPlotter.exe: src/gdjControlPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjControlPlotter.C -o bin/gdjControlPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjResponsePlotter.exe: src/gdjResponsePlotter.C
	$(CXX) $(CXXFLAGS) src/gdjResponsePlotter.C -o bin/gdjResponsePlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjDataMCRawPlotter.exe: src/gdjDataMCRawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjDataMCRawPlotter.C -o bin/gdjDataMCRawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/grlToTex.exe: src/grlToTex.C
	$(CXX) $(CXXFLAGS) src/grlToTex.C -o bin/grlToTex.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/testKeyHandler.exe: src/testKeyHandler.C
	$(CXX) $(CXXFLAGS) src/testKeyHandler.C -o bin/testKeyHandler.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/testSampleHandler.exe: src/testSampleHandler.C
	$(CXX) $(CXXFLAGS) src/testSampleHandler.C -o bin/testSampleHandler.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjToyMultiMix.exe: src/gdjToyMultiMix.C
	$(CXX) $(CXXFLAGS) src/gdjToyMultiMix.C -o bin/gdjToyMultiMix.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotToy.exe: src/gdjPlotToy.C
	$(CXX) $(CXXFLAGS) src/gdjPlotToy.C -o bin/gdjPlotToy.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHistToUnfold.exe: src/gdjHistToUnfold.C
	$(CXX) $(CXXFLAGS) src/gdjHistToUnfold.C -o bin/gdjHistToUnfold.exe $(ROOT) $(INCLUDE) $(LIB) $(ROOUNFOLDLIB) -lATLASGDJ

bin/gdjPlotUnfoldDiagnostics.exe: src/gdjPlotUnfoldDiagnostics.C
	$(CXX) $(CXXFLAGS) src/gdjPlotUnfoldDiagnostics.C -o bin/gdjPlotUnfoldDiagnostics.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotResults.exe: src/gdjPlotResults.C
	$(CXX) $(CXXFLAGS) src/gdjPlotResults.C -o bin/gdjPlotResults.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

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
	rm -f input/ntupleToHist/*~
	rm -f input/histToUnfold/*~
	rm -f input/plotRes/*~
	rm -f lib/*.so
	rm -rf lib
	rm -f obj/*.o
	rm -rf obj
	rm -f src/*~
	rm -f src/#*#
