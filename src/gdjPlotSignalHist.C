//Author: Chris McGinn (2020.08.07)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"

//Local
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/getLinBins.h"
#include "include/globalDebugHandler.h"
#include "include/stringUtil.h"

int gdjPlotSignalHist(std::string inConfigFileName)
{
  globalDebugHandler gBug;
  const bool doGlobalDebug = gBug.GetDoGlobalDebug();

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* inConfig_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> reqParams = {"INFILENAME"};

  if(!checkEnvForParams(inConfig_p, reqParams)) return 1;

  const std::string inFileName = inConfig_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("pdfDir/" + dateStr);

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TEnv* fileConfig_p = (TEnv*)inFile_p->Get("config");
  std::vector<std::string> reqFileParams = {"NGAMMAPTBINS",
					    "GAMMAPTBINSLOW",
					    "GAMMAPTBINSHIGH"};

  if(!checkEnvForParams(fileConfig_p, reqFileParams)) return 1;

  const Int_t nMaxBins = 500;
  Int_t nGammaPtBins = fileConfig_p->GetValue("NGAMMAPTBINS", -1);
  Float_t gammaPtBinsLow = fileConfig_p->GetValue("NGAMMAPTBINSLOW", -1);
  Float_t gammaPtBinsHigh = fileConfig_p->GetValue("NGAMMAPTBINSHIGH", -1);
  Double_t gammaPtBins[nMaxBins+1];
  getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  
  
  inFile_p->Close();
  delete inFile_p;
  delete inConfig_p;
  
  std::cout << "GDJPLOTSIGNALHIST COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotSignalHist.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjPlotSignalHist(argv[1]);
  return retVal;
}
