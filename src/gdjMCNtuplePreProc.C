//Author: Chris McGinn (2020.02.26)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"

//Local                                                                                   
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/globalDebugHandler.h"
#include "include/returnFileList.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"

int gdjMCNtuplePreProc(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;

  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  configParser config(inConfigFileName);

  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  std::string inDirStr = config.GetConfigVal("MCPREPROCDIRNAME");
  if(!check.checkDir(inDirStr)){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  std::string outFileName = config.GetConfigVal("OUTFILENAME");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("gammaJetTree_p", "");

  Float_t pthat_;
  Int_t sampleTag_;
  
  outTree_p->Branch("pthat", &pthat_, "pthat/F");
  outTree_p->Branch("sampleTag", &sampleTag_, "sampleTag/I");
  
  std::vector<std::string> fileList = returnFileList(inDirStr, ".root");
  if(fileList.size() == 0){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' contains no root files. return 1" << std::endl;
    return 1;
  }

  std::map<std::string, std::vector<std::string> > configMap;
  configMap["MCPREPROCDIRNAME"] = {inDirStr};

  //Basic pre-processing for output config
  ULong64_t totalNEntries = 0;
  for(auto const & file : fileList){
    TFile* inFile_p = new TFile(file.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    totalNEntries += inTree_p->GetEntries();
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
    
    configParser tempConfig(inConfig_p);
    std::map<std::string, std::string> tempConfigMap = tempConfig.GetConfigMap();

    for(auto const & val : tempConfigMap){
      bool isFound = false;
      std::vector<std::string> tempVect = configMap[val.first];
      for(unsigned int vI = 0; vI < tempVect.size(); ++vI){
	if(isStrSame(tempVect[vI], val.second)){
	  isFound = true;
	  break;
	}
      }

      if(!isFound){
	tempVect.push_back(val.second);
	configMap[val.first] = tempVect;
      }      
    }
    
    inFile_p->Close();
    delete inFile_p;
  }

  sampleHandler sHandler;
  
  ULong64_t nDiv = TMath::Max((ULong64_t)1, totalNEntries/20);
  ULong64_t currTotalEntries = 0;
  UInt_t nFile = 0;
  for(auto const & file : fileList){
    TFile* inFile_p = new TFile(file.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
    configParser tempConfig(inConfig_p);
    std::string inDataSetName = tempConfig.GetConfigVal("INDATASET");
    if(inDataSetName.size() == 0 || !sHandler.Init(inDataSetName)){
      std::cout << "GDJMCNTUPLEPREPROC ERROR - Given input \'" << file << "\' contains INDATASET \'" << inDataSetName << "\' that is not valid. return 1" << std::endl;

      inFile_p->Close();
      delete inFile_p;

      outFile_p->Close();
      delete outFile_p;
      
      return 1;
    }

    sampleTag_ = sHandler.GetTag();
    
    inTree_p->SetBranchAddress("pthat", &pthat_);

    const ULong64_t nEntries = inTree_p->GetEntries();
    for(ULong64_t entry = 0; entry < nEntries; ++entry){
      if(currTotalEntries%nDiv == 0) std::cout << " Entry " << currTotalEntries << "/" << totalNEntries << "... (File " << nFile << "/" << fileList.size() << ")"  << std::endl;
      inTree_p->GetEntry(entry);

      outTree_p->Fill();
      ++currTotalEntries;
    }

    inFile_p->Close();
    delete inFile_p;

    ++nFile;
  }
  
  outFile_p->cd();

  outTree_p->Write("", TObject::kOverwrite);
  delete outTree_p;
  
  TEnv outConfig;
  for(auto const & val : configMap){  
    for(unsigned int vI = 0; vI < val.second.size(); ++vI){
      std::string configVal = val.first;
      if(val.second.size() != 1) configVal = configVal + "_" + std::to_string(vI);
      outConfig.SetValue(configVal.c_str(), val.second[vI].c_str());
    }    
  }

  outConfig.Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;
  
  std::cout << "GDJMCNTUPLEPREPROC COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjMCNtuplePreProc.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += gdjMCNtuplePreProc(argv[1]);
  return retVal;
}
