//Author: Chris McGinn; Created 2023.02.23
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TTree.h"

//Local
#include "include/checkMakeDir.h"
#include "include/stringUtil.h"

int gdjGetCentMixWeights(const std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  config_p->SetValue("CONFIGFILENAME", inConfigFileName.c_str());

  std::vector<std::string> necessaryParams = {"INFILENAME",
					      "OUTFILENAME"};

  const std::string inFileName = config_p->GetValue("INFILENAME", "");
  if(!check.checkFileExt(inFileName, ".root")) return 1;

  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output"); // check output dir exists; if not create
  check.doCheckMakeDir("output/" + dateStr); // check dated output subdir exists; if not create

  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".txt";
  
  std::map<int, unsigned long long> eventsPerCent, eventsPerCentWith20GeVR4Jet;
  const int centMin = 0;
  const int centMax = 80;
  for(int centIter = centMin; centIter <= centMax; ++centIter){
    eventsPerCent[centIter] = 0;
    eventsPerCentWith20GeVR4Jet[centIter] = 0;
  }

  Int_t cent_;
  std::vector<float>* akt4hi_insitu_jet_pt_p=nullptr;

  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");

  inTree_p->SetBranchStatus("*", 0);
  inTree_p->SetBranchStatus("cent", 1);
  inTree_p->SetBranchStatus("akt4hi_insitu_jet_pt", 1);

  inTree_p->SetBranchAddress("cent", &cent_);
  inTree_p->SetBranchAddress("akt4hi_insitu_jet_pt", &akt4hi_insitu_jet_pt_p);

  const ULong64_t nEntries = inTree_p->GetEntries();
  //Define nDiv 
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..."<< std::endl;
    inTree_p->GetEntry(entry);

    bool oneGoodJet = false;
    for(unsigned int jI = 0; jI < akt4hi_insitu_jet_pt_p->size(); ++jI){
      if(akt4hi_insitu_jet_pt_p->at(jI) > 20.0){
	oneGoodJet = true;
	break;
      }
    }

    ++eventsPerCent[cent_];
    if(oneGoodJet) ++eventsPerCentWith20GeVR4Jet[cent_];
  }
  std::cout << "Processing complete!" << std::endl;

  inFile_p->Close();
  delete inFile_p;

  std::ofstream outFile(outFileName.c_str());
  

  for(auto const& val : eventsPerCent){
    //Commented out initial section which was just checking output values
    /*
    double valueToPrint = eventsPerCentWith20GeVR4Jet[val.first];
    valueToPrint /= (double)val.second;

    std::cout << val.first << ": " << eventsPerCentWith20GeVR4Jet[val.first] << "/" << val.second << "=" << valueToPrint << std::endl;
    */

    outFile << val.first << ", " << eventsPerCentWith20GeVR4Jet[val.first] << ", " << val.second << std::endl;
  }

  outFile.close();

  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjGetCentMixWeights.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjGetCentMixWeights(argv[1]);
  return retVal;
}
