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

HEPMC=/home/cfm/Packages/HepMC2/hepmc-install

ROOUNFOLDDIR=/home/cfm/Packages/RooUnfold/RooUnfold-build/
INCLUDE=-I$(GDJDIR) -I$(ROOUNFOLDDIR) -I$(HEPMC)/include
LIB=-L$(GDJDIR)/lib  -L$(HEPMC)/lib -lHepMC
ROOUNFOLDLIB=-L$(ROOUNFOLDDIR) -lRooUnfold

ROOT=`root-config --cflags --glibs`
FASTJET=`fastjet-config --cxxflags --libs --plugins`

MKDIR_BIN=mkdir -p $(GDJDIR)/bin
MKDIR_LIB=mkdir -p $(GDJDIR)/lib
MKDIR_OBJ=mkdir -p $(GDJDIR)/obj
MKDIR_OUTPUT=mkdir -p $(GDJDIR)/output
MKDIR_PDF=mkdir -p $(GDJDIR)/pdfDir

all: mkdirBin mkdirLib mkdirObj mkdirOutput mkdirPdf obj/binFlattener.o obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o obj/mixMachine.o lib/libATLASGDJ.so bin/gdjNtuplePreProc.exe bin/gdjAnalyzeTxtOut.exe bin/gdjToyMultiMix.exe bin/gdjPlotToy.exe bin/gdjNTupleToHist.exe bin/gdjNTupleToMBHist.exe bin/gdjHistDumper.exe bin/gdjGammaJetResponsePlot.exe bin/gdjMixedEventPlotter.exe bin/gdjPurityPlotter.exe bin/gdjControlPlotter.exe bin/gdjResponsePlotter.exe bin/gdjDataMCRawPlotter.exe bin/gdjPbPbOverPPRawPlotter.exe bin/gdjRCPRawPlotter.exe bin/gdjR4OverR2RawPlotter.exe bin/grlToTex.exe bin/testKeyHandler.exe bin/testSampleHandler.exe bin/gdjPlotMBHist.exe bin/gdjHistToUnfold.exe bin/gdjHistToGenVarPlots.exe bin/gdjPlotUnfoldReweight.exe bin/gdjPlotUnfoldDiagnostics.exe bin/gdjPlotResults.exe bin/gdjHistDQM.exe bin/gdjHEPMCToRoot.exe bin/gdjHEPMCAna.exe bin/gdjHEPMCPlot.exe bin/gdjHEPMCCalib.exe bin/gdjHEPMCCalibPlot.exe bin/gdjRunStabilityPlotter.exe bin/gdjPlotJetVarResponse.exe
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

obj/binFlattener.o: src/binFlattener.C
	$(CXX) $(CXXFLAGS) -fPIC -c src/binFlattener.C -o obj/binFlattener.o $(INCLUDE) $(ROOT)

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
	$(CXX) $(CXXFLAGS) -fPIC -shared -o lib/libATLASGDJ.so obj/binFlattener.o obj/centralityFromInput.o obj/checkMakeDir.o obj/configParser.o obj/globalDebugHandler.o obj/keyHandler.o obj/sampleHandler.o obj/mixMachine.o $(ROOT) $(INCLUDE)

bin/gdjNtuplePreProc.exe: src/gdjNtuplePreProc.C
	$(CXX) $(CXXFLAGS) src/gdjNtuplePreProc.C -o bin/gdjNtuplePreProc.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjAnalyzeTxtOut.exe: src/gdjAnalyzeTxtOut.C
	$(CXX) $(CXXFLAGS) src/gdjAnalyzeTxtOut.C -o bin/gdjAnalyzeTxtOut.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

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

bin/gdjPurityPlotter.exe: src/gdjPurityPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjPurityPlotter.C -o bin/gdjPurityPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjControlPlotter.exe: src/gdjControlPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjControlPlotter.C -o bin/gdjControlPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjResponsePlotter.exe: src/gdjResponsePlotter.C
	$(CXX) $(CXXFLAGS) src/gdjResponsePlotter.C -o bin/gdjResponsePlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjDataMCRawPlotter.exe: src/gdjDataMCRawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjDataMCRawPlotter.C -o bin/gdjDataMCRawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPbPbOverPPRawPlotter.exe: src/gdjPbPbOverPPRawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjPbPbOverPPRawPlotter.C -o bin/gdjPbPbOverPPRawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjRCPRawPlotter.exe: src/gdjRCPRawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjRCPRawPlotter.C -o bin/gdjRCPRawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjR4OverR2RawPlotter.exe: src/gdjR4OverR2RawPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjR4OverR2RawPlotter.C -o bin/gdjR4OverR2RawPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

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

bin/gdjHistToGenVarPlots.exe: src/gdjHistToGenVarPlots.C
	$(CXX) $(CXXFLAGS) src/gdjHistToGenVarPlots.C -o bin/gdjHistToGenVarPlots.exe $(ROOT) $(INCLUDE) $(LIB) $(ROOUNFOLDLIB) -lATLASGDJ

bin/gdjPlotUnfoldReweight.exe: src/gdjPlotUnfoldReweight.C
	$(CXX) $(CXXFLAGS) src/gdjPlotUnfoldReweight.C -o bin/gdjPlotUnfoldReweight.exe $(ROOT) $(INCLUDE) $(LIB) $(ROOUNFOLDLIB) -lATLASGDJ

bin/gdjPlotUnfoldDiagnostics.exe: src/gdjPlotUnfoldDiagnostics.C
	$(CXX) $(CXXFLAGS) src/gdjPlotUnfoldDiagnostics.C -o bin/gdjPlotUnfoldDiagnostics.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotResults.exe: src/gdjPlotResults.C
	$(CXX) $(CXXFLAGS) src/gdjPlotResults.C -o bin/gdjPlotResults.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHistDQM.exe: src/gdjHistDQM.C
	$(CXX) $(CXXFLAGS) src/gdjHistDQM.C -o bin/gdjHistDQM.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHEPMCToRoot.exe: src/gdjHEPMCToRoot.C
	$(CXX) $(CXXFLAGS) src/gdjHEPMCToRoot.C -o bin/gdjHEPMCToRoot.exe $(ROOT) $(FASTJET) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHEPMCAna.exe: src/gdjHEPMCAna.C
	$(CXX) $(CXXFLAGS) src/gdjHEPMCAna.C -o bin/gdjHEPMCAna.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHEPMCCalib.exe: src/gdjHEPMCCalib.C
	$(CXX) $(CXXFLAGS) src/gdjHEPMCCalib.C -o bin/gdjHEPMCCalib.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHEPMCPlot.exe: src/gdjHEPMCPlot.C
	$(CXX) $(CXXFLAGS) src/gdjHEPMCPlot.C -o bin/gdjHEPMCPlot.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjHEPMCCalibPlot.exe: src/gdjHEPMCCalibPlot.C
	$(CXX) $(CXXFLAGS) src/gdjHEPMCCalibPlot.C -o bin/gdjHEPMCCalibPlot.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjRunStabilityPlotter.exe: src/gdjRunStabilityPlotter.C
	$(CXX) $(CXXFLAGS) src/gdjRunStabilityPlotter.C -o bin/gdjRunStabilityPlotter.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

bin/gdjPlotJetVarResponse.exe: src/gdjPlotJetVarResponse.C
	$(CXX) $(CXXFLAGS) src/gdjPlotJetVarResponse.C -o bin/gdjPlotJetVarResponse.exe $(ROOT) $(INCLUDE) $(LIB) -lATLASGDJ

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
