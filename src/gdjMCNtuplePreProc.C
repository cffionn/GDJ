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
#include "include/centralityFromInput.h"
#include "include/checkMakeDir.h"
#include "include/configParser.h"
#include "include/globalDebugHandler.h"
#include "include/returnFileList.h"
#include "include/sampleHandler.h"
#include "include/stringUtil.h"
#include "include/treeUtil.h"

int gdjMCNtuplePreProc(std::string inConfigFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".txt")) return 1;
  
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  configParser config(inConfigFileName);

  std::vector<std::string> necessaryParams = {"MCPREPROCDIRNAME",
					      "OUTFILENAME",
					      "CENTFILENAME",
					      "ISPP"};

  if(!config.ContainsParamSet(necessaryParams)) return 1;
  
  if(doGlobalDebug) std::cout << "GLOBAL DEBUG FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

  const std::string inDirStr = config.GetConfigVal("MCPREPROCDIRNAME");
  const std::string inCentFileName = config.GetConfigVal("CENTFILENAME");
  if(!check.checkDir(inDirStr)){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' is not valid. return 1" << std::endl;
    return 1;
  }
  if(!check.checkFileExt(inCentFileName, "txt")) return 1;

  centralityFromInput centTable(inCentFileName);
  if(doGlobalDebug) centTable.PrintTableTex();
  
  const std::string dateStr = getDateStr();
  check.doCheckMakeDir("output");
  check.doCheckMakeDir("output/" + dateStr);

  std::string outFileName = config.GetConfigVal("OUTFILENAME");
  if(outFileName.find(".") != std::string::npos) outFileName = outFileName.substr(0, outFileName.rfind("."));
  outFileName = "output/" + dateStr + "/" + outFileName + "_" + dateStr + ".root";

  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  TTree* outTree_p = new TTree("gammaJetTree_p", "");

  Int_t runNumber_, eventNumber_;
  UInt_t lumiBlock_;
  Bool_t passesToroid_;

  Float_t pthat_;
  Int_t cent_;
  Int_t sampleTag_;
  Float_t sampleWeight_;
  Float_t ncollWeight_;
  Float_t fullWeight_;
  
  const Int_t nTreeParton_ = 2;
  Float_t treePartonPt_[nTreeParton_];
  Float_t treePartonEta_[nTreeParton_];
  Float_t treePartonPhi_[nTreeParton_];
  Int_t treePartonId_[nTreeParton_];

  Float_t actualInteractionsPerCrossing_, averageInteractionsPerCrossing_;
  Int_t nvert_;
  std::vector<float>* vert_x_p=nullptr;
  std::vector<float>* vert_y_p=nullptr;
  std::vector<float>* vert_z_p=nullptr;
  std::vector<int>* vert_type_p=nullptr;
  std::vector<int>* vert_ntrk_p=nullptr;

  Float_t fcalA_et_, fcalC_et_;

  Int_t truth_n_;
  std::vector<float>* truth_charge_p=nullptr;
  std::vector<float>* truth_pt_p=nullptr;
  std::vector<float>* truth_eta_p=nullptr;
  std::vector<float>* truth_phi_p=nullptr;

  Int_t akt2hi_jet_n_;
  std::vector<float>* akt2hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt2hi_constit_xcalib_jet_e_p=nullptr;
  std::vector<int>* akt2hi_truthpos_p=nullptr;
  
  Int_t akt4hi_jet_n_;
  std::vector<float>* akt4hi_em_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_em_xcalib_jet_e_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_pt_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_eta_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_phi_p=nullptr;
  std::vector<float>* akt4hi_constit_xcalib_jet_e_p=nullptr;
  std::vector<int>* akt4hi_truthpos_p=nullptr;

  Int_t photon_n_;
  std::vector<float>* photon_pt_p=nullptr;
  std::vector<float>* photon_eta_p=nullptr;
  std::vector<float>* photon_phi_p=nullptr;
  std::vector<bool>* photon_tight_p=nullptr;
  std::vector<bool>* photon_loose_p=nullptr;
  std::vector<unsigned int>* photon_isem_p=nullptr;
  std::vector<int>* photon_convFlag_p=nullptr;
  std::vector<float>* photon_Rconv_p=nullptr;
  std::vector<float>* photon_etcone20_p=nullptr;
  std::vector<float>* photon_etcone30_p=nullptr;
  std::vector<float>* photon_etcone40_p=nullptr;
  std::vector<float>* photon_topoetcone20_p=nullptr;
  std::vector<float>* photon_topoetcone30_p=nullptr;
  std::vector<float>* photon_topoetcone40_p=nullptr;
  std::vector<float>* photon_Rhad1_p=nullptr;
  std::vector<float>* photon_Rhad_p=nullptr;
  std::vector<float>* photon_e277_p=nullptr;
  std::vector<float>* photon_Reta_p=nullptr;
  std::vector<float>* photon_Rphi_p=nullptr;
  std::vector<float>* photon_weta1_p=nullptr;
  std::vector<float>* photon_weta2_p=nullptr;
  std::vector<float>* photon_wtots1_p=nullptr;
  std::vector<float>* photon_f1_p=nullptr;
  std::vector<float>* photon_f3_p=nullptr;
  std::vector<float>* photon_fracs1_p=nullptr;
  std::vector<float>* photon_DeltaE_p=nullptr;
  std::vector<float>* photon_Eratio_p=nullptr;
  Int_t akt2_truth_jet_n_;

  std::vector<float>* akt2_truth_jet_pt_p=nullptr;
  std::vector<float>* akt2_truth_jet_eta_p=nullptr;
  std::vector<float>* akt2_truth_jet_phi_p=nullptr;
  std::vector<float>* akt2_truth_jet_e_p=nullptr;
  std::vector<int>* akt2_truth_jet_recopos_p=nullptr;
  Int_t akt4_truth_jet_n_;
  std::vector<float>* akt4_truth_jet_pt_p=nullptr;
  std::vector<float>* akt4_truth_jet_eta_p=nullptr;
  std::vector<float>* akt4_truth_jet_phi_p=nullptr;
  std::vector<float>* akt4_truth_jet_e_p=nullptr;
  std::vector<int>* akt4_truth_jet_recopos_p=nullptr;
  
  outTree_p->Branch("runNumber", &runNumber_, "runNumber/I");
  outTree_p->Branch("eventNumber", &eventNumber_, "eventNumber/I");
  outTree_p->Branch("lumiBlock", &lumiBlock_, "lumiBlock/i");
  outTree_p->Branch("passesToroid", &passesToroid_, "passesToroid/i");

  outTree_p->Branch("pthat", &pthat_, "pthat/F");
  outTree_p->Branch("cent", &cent_, "cent/I");
  outTree_p->Branch("sampleTag", &sampleTag_, "sampleTag/I");
  outTree_p->Branch("sampleWeight", &sampleWeight_, "sampleWeight/F");
  outTree_p->Branch("ncollWeight", &ncollWeight_, "ncollWeight/F");
  outTree_p->Branch("fullWeight", &fullWeight_, "fullWeight/F");

  outTree_p->Branch("treePartonPt", treePartonPt_, ("treePartonPt[" + std::to_string(nTreeParton_) + "]/F").c_str());
  outTree_p->Branch("treePartonEta", treePartonEta_, ("treePartonEta[" + std::to_string(nTreeParton_) + "]/F").c_str());
  outTree_p->Branch("treePartonPhi", treePartonPhi_, ("treePartonPhi[" + std::to_string(nTreeParton_) + "]/F").c_str());
  outTree_p->Branch("treePartonId", treePartonId_, ("treePartonId[" + std::to_string(nTreeParton_) + "]/I").c_str());
  
  outTree_p->Branch("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_, "actualInteractionsPerCrossing/F");
  outTree_p->Branch("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_, "averageInteractionsPerCrossing/F");
  outTree_p->Branch("nvert", &nvert_, "nvert/I");

  outTree_p->Branch("vert_x", &vert_x_p);
  outTree_p->Branch("vert_y", &vert_y_p);
  outTree_p->Branch("vert_z", &vert_z_p);
  outTree_p->Branch("vert_type", &vert_type_p);
  outTree_p->Branch("vert_ntrk", &vert_ntrk_p);

  outTree_p->Branch("fcalA_et", &fcalA_et_, "fcalA_et/F");
  outTree_p->Branch("fcalC_et", &fcalC_et_, "fcalC_et/F");

  outTree_p->Branch("truth_n", &truth_n_, "truth_n/I");
  outTree_p->Branch("truth_charge", &truth_charge_p);
  outTree_p->Branch("truth_pt", &truth_pt_p);
  outTree_p->Branch("truth_eta", &truth_eta_p);
  outTree_p->Branch("truth_phi", &truth_phi_p);

  outTree_p->Branch("akt2hi_jet_n", &akt2hi_jet_n_, "akt2hi_jet_n/I");
  outTree_p->Branch("akt2hi_em_xcalib_jet_pt", &akt2hi_em_xcalib_jet_pt_p);
  outTree_p->Branch("akt2hi_em_xcalib_jet_eta", &akt2hi_em_xcalib_jet_eta_p);
  outTree_p->Branch("akt2hi_em_xcalib_jet_phi", &akt2hi_em_xcalib_jet_phi_p);
  outTree_p->Branch("akt2hi_em_xcalib_jet_e", &akt2hi_em_xcalib_jet_e_p);
  outTree_p->Branch("akt2hi_constit_xcalib_jet_pt", &akt2hi_constit_xcalib_jet_pt_p);
  outTree_p->Branch("akt2hi_constit_xcalib_jet_eta", &akt2hi_constit_xcalib_jet_eta_p);
  outTree_p->Branch("akt2hi_constit_xcalib_jet_phi", &akt2hi_constit_xcalib_jet_phi_p);
  outTree_p->Branch("akt2hi_constit_xcalib_jet_e", &akt2hi_constit_xcalib_jet_e_p);
  outTree_p->Branch("akt2hi_truthpos", &akt2hi_truthpos_p);

  outTree_p->Branch("akt4hi_jet_n", &akt4hi_jet_n_, "akt4hi_jet_n/I");
  outTree_p->Branch("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
  outTree_p->Branch("akt4hi_em_xcalib_jet_e", &akt4hi_em_xcalib_jet_e_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_pt", &akt4hi_constit_xcalib_jet_pt_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_eta", &akt4hi_constit_xcalib_jet_eta_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_phi", &akt4hi_constit_xcalib_jet_phi_p);
  outTree_p->Branch("akt4hi_constit_xcalib_jet_e", &akt4hi_constit_xcalib_jet_e_p);
  outTree_p->Branch("akt4hi_truthpos", &akt4hi_truthpos_p);

  outTree_p->Branch("photon_n", &photon_n_, "photon_n/I");
  outTree_p->Branch("photon_pt", &photon_pt_p);
  outTree_p->Branch("photon_eta", &photon_eta_p);
  outTree_p->Branch("photon_phi", &photon_phi_p);
  outTree_p->Branch("photon_tight", &photon_tight_p);
  outTree_p->Branch("photon_loose", &photon_loose_p);
  outTree_p->Branch("photon_isem", &photon_isem_p);
  outTree_p->Branch("photon_convFlag", &photon_convFlag_p);
  outTree_p->Branch("photon_Rconv", &photon_Rconv_p);
  outTree_p->Branch("photon_etcone20", &photon_etcone20_p);
  outTree_p->Branch("photon_etcone30", &photon_etcone30_p);
  outTree_p->Branch("photon_etcone40", &photon_etcone40_p);
  outTree_p->Branch("photon_topoetcone20", &photon_topoetcone20_p);
  outTree_p->Branch("photon_topoetcone30", &photon_topoetcone30_p);
  outTree_p->Branch("photon_topoetcone40", &photon_topoetcone40_p);
  outTree_p->Branch("photon_Rhad1", &photon_Rhad1_p);
  outTree_p->Branch("photon_Rhad", &photon_Rhad_p);
  outTree_p->Branch("photon_e277", &photon_e277_p);
  outTree_p->Branch("photon_Reta", &photon_Reta_p);
  outTree_p->Branch("photon_Rphi", &photon_Rphi_p);
  outTree_p->Branch("photon_weta1", &photon_weta1_p);
  outTree_p->Branch("photon_weta2", &photon_weta2_p);
  outTree_p->Branch("photon_weta2", &photon_weta2_p);
  outTree_p->Branch("photon_wtots1", &photon_wtots1_p);
  outTree_p->Branch("photon_f1", &photon_f1_p);
  outTree_p->Branch("photon_f3", &photon_f3_p);
  outTree_p->Branch("photon_fracs1", &photon_fracs1_p);
  outTree_p->Branch("photon_DeltaE", &photon_DeltaE_p);
  outTree_p->Branch("photon_Eratio", &photon_Eratio_p);

  outTree_p->Branch("akt2_truth_jet_n", &akt2_truth_jet_n_, "akt2_truth_jet_n/I");
  outTree_p->Branch("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
  outTree_p->Branch("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
  outTree_p->Branch("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
  outTree_p->Branch("akt2_truth_jet_e", &akt2_truth_jet_e_p);
  outTree_p->Branch("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);

  outTree_p->Branch("akt4_truth_jet_n", &akt4_truth_jet_n_, "akt4_truth_jet_n/I");
  outTree_p->Branch("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
  outTree_p->Branch("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
  outTree_p->Branch("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
  outTree_p->Branch("akt4_truth_jet_e", &akt4_truth_jet_e_p);
  outTree_p->Branch("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);

  std::vector<std::string> fileList = returnFileList(inDirStr, ".root");
  if(fileList.size() == 0){
    std::cout << "GDJMCNTUPLEPREPROC ERROR - Given MCPREPROCDIRNAME \'" << inDirStr << "\' in config \'" << inConfigFileName << "\' contains no root files. return 1" << std::endl;
    return 1;
  }

  std::map<int, unsigned long long> tagToCounts;
  std::map<int, double> tagToXSec;
  std::map<int, double> tagToFilterEff;
  std::vector<int> minPthats;
  std::map<int, int> minPthatToTag;
  std::map<int, double> tagToWeight;
  
  std::map<std::string, std::vector<std::string> > configMap;
  configMap["MCPREPROCDIRNAME"] = {inDirStr};

  //Basic pre-processing for output config
  std::vector<std::string> listOfBranchesIn;
  std::vector<std::string> listOfBranchesOut = getVectBranchList(outTree_p);
  ULong64_t totalNEntries = 0;
  for(auto const & file : fileList){
    TFile* inFile_p = new TFile(file.c_str(), "READ");
    TTree* inTree_p = (TTree*)inFile_p->Get("gammaJetTree_p");
    totalNEntries += inTree_p->GetEntries();
    TEnv* inConfig_p = (TEnv*)inFile_p->Get("config");
    
    configParser tempConfig(inConfig_p);
    std::map<std::string, std::string> tempConfigMap = tempConfig.GetConfigMap();

    sampleHandler sHandler;
    
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
    
    if(tagToCounts.count(sampleTag_) == 0){
      tagToCounts[sampleTag_] = inTree_p->GetEntries();
      tagToXSec[sampleTag_] = sHandler.GetXSection();
      tagToFilterEff[sampleTag_] = sHandler.GetFilterEff();
      int minPthat = sHandler.GetMinPthat();
      minPthats.push_back(minPthat);
      minPthatToTag[sHandler.GetMinPthat()] = sampleTag_;
    }
    else tagToCounts[sampleTag_] += inTree_p->GetEntries();
    
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

    std::vector<std::string> tempBranches = getVectBranchList(inTree_p);
    if(listOfBranchesIn.size() == 0) listOfBranchesIn = tempBranches;
    else{
      for(auto const & branch : tempBranches){
	if(!vectContainsStr(branch, &listOfBranchesIn)){
	  listOfBranchesIn.push_back(branch);
	}
      }
    }
        
    inFile_p->Close();
    delete inFile_p;
  }

  std::sort(std::begin(minPthats), std::end(minPthats));
  std::cout << "Min pthats: " << std::endl;
  for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
    int tag = minPthatToTag[minPthats[pI]];
    double weight = tagToXSec[tag]*tagToFilterEff[tag]/(double)tagToCounts[tag];
    tagToWeight[tag] = weight;
  }

  for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
    if(pI == 0) continue;

    int tag = minPthatToTag[minPthats[pI]];
    int tag0 = minPthatToTag[minPthats[0]];

    tagToWeight[tag] /= tagToWeight[tag0];
  }
  
  int tag0 = minPthatToTag[minPthats[0]];
  tagToWeight[tag0] = 1.0;
  
  for(unsigned int pI = 0; pI < minPthats.size(); ++pI){
    int tag = minPthatToTag[minPthats[pI]];

    std::cout << " " << pI << "/" << minPthats.size() << ": " << minPthats[pI] << ", " << tagToWeight[tag] << std::endl;
  }

  bool allBranchesGood = true;
  for(auto const & branchIn : listOfBranchesIn){
    bool containsBranch = vectContainsStr(branchIn, &listOfBranchesOut);
    allBranchesGood = allBranchesGood && containsBranch;

    if(!containsBranch) std::cout << "GDJMCNTUPLEPREPROC ERROR - \'" << branchIn << "\' branch is not part of output tree. please fix macro, return 1" << std::endl;
  }

  if(!allBranchesGood){
    outFile_p->Close();
    delete outFile_p;
    return 1;
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
    sampleWeight_ = tagToWeight[sampleTag_];
   
    inTree_p->SetBranchAddress("runNumber", &runNumber_);
    inTree_p->SetBranchAddress("eventNumber", &eventNumber_);
    inTree_p->SetBranchAddress("lumiBlock", &lumiBlock_);
    inTree_p->SetBranchAddress("pthat", &pthat_);

    inTree_p->SetBranchAddress("treePartonPt", treePartonPt_);
    inTree_p->SetBranchAddress("treePartonEta", treePartonEta_);
    inTree_p->SetBranchAddress("treePartonPhi", treePartonPhi_);
    inTree_p->SetBranchAddress("treePartonId", treePartonId_);

    inTree_p->SetBranchAddress("actualInteractionsPerCrossing", &actualInteractionsPerCrossing_);
    inTree_p->SetBranchAddress("averageInteractionsPerCrossing", &averageInteractionsPerCrossing_);

    inTree_p->SetBranchAddress("nvert", &nvert_);
    inTree_p->SetBranchAddress("vert_x", &vert_x_p);
    inTree_p->SetBranchAddress("vert_y", &vert_y_p);
    inTree_p->SetBranchAddress("vert_z", &vert_z_p);
    inTree_p->SetBranchAddress("vert_type", &vert_type_p);
    inTree_p->SetBranchAddress("vert_ntrk", &vert_ntrk_p);

    inTree_p->SetBranchAddress("fcalA_et", &fcalA_et_);
    inTree_p->SetBranchAddress("fcalC_et", &fcalC_et_);

    inTree_p->SetBranchAddress("truth_n", &truth_n_);
    inTree_p->SetBranchAddress("truth_charge", &truth_charge_p);
    inTree_p->SetBranchAddress("truth_pt", &truth_pt_p);
    inTree_p->SetBranchAddress("truth_eta", &truth_eta_p);
    inTree_p->SetBranchAddress("truth_phi", &truth_phi_p);

    inTree_p->SetBranchAddress("akt2hi_jet_n", &akt2hi_jet_n_);
    inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_pt", &akt2hi_em_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_eta", &akt2hi_em_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_phi", &akt2hi_em_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt2hi_em_xcalib_jet_e", &akt2hi_em_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_pt", &akt2hi_constit_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_eta", &akt2hi_constit_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_phi", &akt2hi_constit_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt2hi_constit_xcalib_jet_e", &akt2hi_constit_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt2hi_truthpos", &akt2hi_truthpos_p);

    inTree_p->SetBranchAddress("akt4hi_jet_n", &akt4hi_jet_n_);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_pt", &akt4hi_em_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_eta", &akt4hi_em_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_phi", &akt4hi_em_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_em_xcalib_jet_e", &akt4hi_em_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_pt", &akt4hi_constit_xcalib_jet_pt_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_eta", &akt4hi_constit_xcalib_jet_eta_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_phi", &akt4hi_constit_xcalib_jet_phi_p);
    inTree_p->SetBranchAddress("akt4hi_constit_xcalib_jet_e", &akt4hi_constit_xcalib_jet_e_p);
    inTree_p->SetBranchAddress("akt4hi_truthpos", &akt4hi_truthpos_p);

    inTree_p->SetBranchAddress("photon_n", &photon_n_);
    inTree_p->SetBranchAddress("photon_pt", &photon_pt_p);
    inTree_p->SetBranchAddress("photon_eta", &photon_eta_p);
    inTree_p->SetBranchAddress("photon_phi", &photon_phi_p);
    inTree_p->SetBranchAddress("photon_tight", &photon_tight_p);
    inTree_p->SetBranchAddress("photon_loose", &photon_loose_p);
    inTree_p->SetBranchAddress("photon_isem", &photon_isem_p);
    inTree_p->SetBranchAddress("photon_convFlag", &photon_convFlag_p);
    inTree_p->SetBranchAddress("photon_Rconv", &photon_Rconv_p);
    inTree_p->SetBranchAddress("photon_etcone20", &photon_etcone20_p);
    inTree_p->SetBranchAddress("photon_etcone30", &photon_etcone30_p);
    inTree_p->SetBranchAddress("photon_etcone40", &photon_etcone40_p);
    inTree_p->SetBranchAddress("photon_topoetcone20", &photon_topoetcone20_p);
    inTree_p->SetBranchAddress("photon_topoetcone30", &photon_topoetcone30_p);
    inTree_p->SetBranchAddress("photon_topoetcone40", &photon_topoetcone40_p);
    inTree_p->SetBranchAddress("photon_Rhad1", &photon_Rhad1_p);
    inTree_p->SetBranchAddress("photon_Rhad", &photon_Rhad_p);
    inTree_p->SetBranchAddress("photon_e277", &photon_e277_p);
    inTree_p->SetBranchAddress("photon_Reta", &photon_Reta_p);
    inTree_p->SetBranchAddress("photon_Rphi", &photon_Rphi_p);
    inTree_p->SetBranchAddress("photon_weta1", &photon_weta1_p);
    inTree_p->SetBranchAddress("photon_weta2", &photon_weta2_p);
    inTree_p->SetBranchAddress("photon_wtots1", &photon_wtots1_p);
    inTree_p->SetBranchAddress("photon_f1", &photon_f1_p);
    inTree_p->SetBranchAddress("photon_f3", &photon_f3_p);
    inTree_p->SetBranchAddress("photon_fracs1", &photon_fracs1_p);
    inTree_p->SetBranchAddress("photon_DeltaE", &photon_DeltaE_p);
    inTree_p->SetBranchAddress("photon_Eratio", &photon_Eratio_p);

    inTree_p->SetBranchAddress("akt2_truth_jet_n", &akt2_truth_jet_n_);
    inTree_p->SetBranchAddress("akt2_truth_jet_pt", &akt2_truth_jet_pt_p);
    inTree_p->SetBranchAddress("akt2_truth_jet_eta", &akt2_truth_jet_eta_p);
    inTree_p->SetBranchAddress("akt2_truth_jet_phi", &akt2_truth_jet_phi_p);
    inTree_p->SetBranchAddress("akt2_truth_jet_e", &akt2_truth_jet_e_p);
    inTree_p->SetBranchAddress("akt2_truth_jet_recopos", &akt2_truth_jet_recopos_p);

    inTree_p->SetBranchAddress("akt4_truth_jet_n", &akt4_truth_jet_n_);
    inTree_p->SetBranchAddress("akt4_truth_jet_pt", &akt4_truth_jet_pt_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_eta", &akt4_truth_jet_eta_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_phi", &akt4_truth_jet_phi_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_e", &akt4_truth_jet_e_p);
    inTree_p->SetBranchAddress("akt4_truth_jet_recopos", &akt4_truth_jet_recopos_p);

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
