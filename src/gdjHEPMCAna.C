//Author: Chris McGinn (2022.11.21)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

//Local
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/etaPhiFunc.h"
#include "include/getLinBins.h"
#include "include/getLogBins.h"
#include "include/ghostUtil.h"
#include "include/globalDebugHandler.h"
#include "include/histDefUtility.h"
#include "include/photonUtil.h"
#include "include/stringUtil.h"
#include "include/varUtil.h"

int gdjHEPMCAna(std::string inConfigFileName)
{
  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Global debug handler - do 'export DOGLOBALDEBUGROOT=1' at command line to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  
  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  //Get the configuration file defining this job
  TEnv* config_p = new TEnv(inConfigFileName.c_str());

  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INFILENAME",
					      "OUTFILENAME",
					      "INUNFOLDFILENAME"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  //Get params in c++ variables we can use
  std::string inFileName = config_p->GetValue("INFILENAME", "");
  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  
  if(!check.checkFileExt(inFileName, ".root")) return 1;
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;

  //We need to extract some info from the file we intend to eventually compare with jewel
  TFile* unfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");
  TEnv* unfoldConfig_p = (TEnv*)unfoldFile_p->Get("config");
  std::vector<std::string> necessaryUnfoldParams = {"NGAMMAPTBINS",
                                                    "GAMMAPTBINSLOW",
                                                    "GAMMAPTBINSHIGH",
                                                    "GAMMAPTBINSDOLOG",
                                                    "GAMMAPTBINSDOCUSTOM",
						    "JTPTBINSLOWRECO",
						    "JTPTBINSHIGHRECO",
						    "JTETABINSLOW",
						    "JTETABINSHIGH",
						    "GAMMAEXCLUSIONDR",
						    "MIXJETEXCLUSIONDR",
						    "GAMMAJTDPHI",
						    "GAMMAMULTIJTDPHI",
						    "VARNAME"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(unfoldConfig_p, necessaryUnfoldParams)) return 1;

  //Get gammaptbins, after defining a bin max
  const Int_t nMaxPtBins = 200;
  const Int_t nGammaPtBins = unfoldConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = unfoldConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = unfoldConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Float_t gammaPtBinsLowReco = unfoldConfig_p->GetValue("GAMMAPTBINSLOWRECO", 80.0);
  const Float_t gammaPtBinsHighReco = unfoldConfig_p->GetValue("GAMMAPTBINSHIGHRECO", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)unfoldConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)unfoldConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxPtBins+1];

  //We need to propagate the gamma pt bins forward  
  config_p->SetValue("NGAMMAPTBINS", nGammaPtBins);
  config_p->SetValue("GAMMAPTBINSLOW", gammaPtBinsLow);
  config_p->SetValue("GAMMAPTBINSHIGH", gammaPtBinsHigh);
  config_p->SetValue("GAMMAPTBINSLOWRECO", gammaPtBinsLowReco);
  config_p->SetValue("GAMMAPTBINSHIGHRECO", gammaPtBinsHighReco);
  config_p->SetValue("GAMMAPTBINSDOLOG", gammaPtBinsDoLog);
  config_p->SetValue("GAMMAPTBINSDOCUSTOM", gammaPtBinsDoCustom);
  config_p->SetValue("GAMMAPTBINSCUSTOM", unfoldConfig_p->GetValue("GAMMAPTBINSCUSTOM", ""));
  
  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;

    std::string gammaPtBinsStr = unfoldConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBins != nGammaPtBinsTemp){
      std::cout << "nGammaPtBins dont match - by stored value \'" << nGammaPtBins << "\', via vector \'\
" << nGammaPtBinsTemp << "\'" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }

    if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  //Need a vector version of the bins for gamma pt bin location
  std::vector<float> gammaPtBinsVect;
  for(Int_t gI = 0; gI < nGammaPtBins+1; ++gI){
    gammaPtBinsVect.push_back(gammaPtBins[gI]);
  }
  
  //Get Reco PT bounds
  const Float_t jtPtBinsLowReco = unfoldConfig_p->GetValue("JTPTBINSLOWRECO", -999.0);
  const Float_t jtPtBinsHighReco = unfoldConfig_p->GetValue("JTPTBINSHIGHRECO", -999.0);
  const Float_t jtEtaBinsLow = unfoldConfig_p->GetValue("JTETABINSLOW", -999.0);
  const Float_t jtEtaBinsHigh = unfoldConfig_p->GetValue("JTETABINSHIGH", -999.0);

  //Get Var name and derivatives
  std::string varName = unfoldConfig_p->GetValue("VARNAME", "");
  std::string varNameUpper = strLowerToUpper(varName);
  std::string varNameLower = returnAllLowercaseString(varName);
  std::string varNameLabel = varNameToLabel(varName);

  //get specific gamma-jet cut values
  const Double_t gammaJtDRExclusionCut = unfoldConfig_p->GetValue("GAMMAEXCLUSIONDR", 100.0);
  const Double_t mixJtDRExclusionCut = unfoldConfig_p->GetValue("MIXJETEXCLUSIONDR", 100.0);
  const Double_t gammaJtDPhiCut = mathStringToNum(unfoldConfig_p->GetValue("GAMMAJTDPHI", ""), doGlobalDebug);
  const Double_t gammaMultiJtDPhiCut = mathStringToNum(unfoldConfig_p->GetValue("GAMMAMULTIJTDPHI", ""), doGlobalDebug);

  //Some more of the above needs to be propagated forward via config
  config_p->SetValue("JTPTBINSLOWRECO", jtPtBinsLowReco);
  config_p->SetValue("JTPTBINSHIGHRECO", jtPtBinsHighReco);
  config_p->SetValue("JTETABINSLOW", jtEtaBinsLow);
  config_p->SetValue("JTETABINSHIGH", jtEtaBinsHigh);

  config_p->SetValue("VARNAME", varName.c_str());
  
  config_p->SetValue("GAMMAEXCLUSIONDR", gammaJtDRExclusionCut);
  config_p->SetValue("MIXJETEXCLUSIONDR", mixJtDRExclusionCut);
  config_p->SetValue("GAMMAJTDPHI", unfoldConfig_p->GetValue("GAMMAJTDPHI", ""));
  config_p->SetValue("GAMMAMULTIJTDPHI", unfoldConfig_p->GetValue("GAMMAMULTIJTDPHI", ""));
  
  //multijet defined by var name
  const Bool_t isMultijet = isStrSame(varNameLower, "xjj") || isStrSame(varNameLower, "dphijjg") || isStrSame(varNameLower, "ajj") || isStrSame(varNameLower, "dphijj") || isStrSame(varNameLower, "drjj");
  

  std::string varPrefix = "";
  if(isStrSame(varNameUpper, "PT")) varPrefix = "JT";
  
  std::string binVarStr = varPrefix + varNameUpper;
  if(isStrSame(binVarStr, "AJJ")) binVarStr = "AJ";
  else if(isStrSame(binVarStr, "DRJJ")) binVarStr = "DR";


  //get var bins
  Int_t nVarBins = unfoldConfig_p->GetValue(("N" + binVarStr + "BINS").c_str(), -1);
  Float_t varBinsLow = unfoldConfig_p->GetValue((binVarStr + "BINSLOW").c_str(), -1.0);
  Float_t varBinsHigh = unfoldConfig_p->GetValue((binVarStr + "BINSHIGH").c_str(), -1.0);
  Bool_t varBinsDoLog = unfoldConfig_p->GetValue((binVarStr + "BINSDOLOG").c_str(), 0);
  Bool_t varBinsDoCustom = unfoldConfig_p->GetValue((binVarStr + "BINSDOCUSTOM").c_str(), 0);
  Float_t varBinsLowReco = unfoldConfig_p->GetValue((binVarStr + "BINSLOWRECO").c_str(), varBinsLow);
  Float_t varBinsHighReco = unfoldConfig_p->GetValue((binVarStr + "BINSHIGHRECO").c_str(), varBinsHigh);
  Double_t varBins[nMaxPtBins+1];

  if(varBinsDoLog) getLogBins(varBinsLow, varBinsHigh, nVarBins, varBins);
  else if(varBinsDoCustom){
    std::string varBinsStr = unfoldConfig_p->GetValue((binVarStr + "BINSCUSTOM").c_str(), "");
    if(varBinsStr.size() == 0){
      std::cout << "No custom bins given. return 1" << std::endl;
      return 1;
    }
    std::vector<float> varBinsTemp = strToVectF(varBinsStr);
    Int_t nVarBinsTemp = ((Int_t)varBinsTemp.size()) - 1;
    if(nVarBins != nVarBinsTemp){
      return 1;
    }

    for(unsigned int gI = 0; gI < varBinsTemp.size(); ++gI){
      varBins[gI] = varBinsTemp[gI];
    }
  }
  else getLinBins(varBinsLow, varBinsHigh, nVarBins, varBins);

  //Need to propagate the varbins as well
  config_p->SetValue(("N" + binVarStr + "BINS").c_str(), nVarBins);
  config_p->SetValue((binVarStr + "BINSLOW").c_str(), varBinsLow);
  config_p->SetValue((binVarStr + "BINSHIGH").c_str(), varBinsHigh);
  config_p->SetValue((binVarStr + "BINSDOLOG").c_str(), varBinsDoLog);
  config_p->SetValue((binVarStr + "BINSDOCUSTOM").c_str(), varBinsDoCustom);
  config_p->SetValue((binVarStr + "BINSCUSTOM").c_str(), unfoldConfig_p->GetValue((binVarStr + "BINSCUSTOM").c_str(), ""));
  
  TFile* inFile_p = new TFile(inFileName.c_str(), "READ");
  TTree* jewelTree_p = (TTree*)inFile_p->Get("jewelTree");
  TEnv* inFileConfig_p = (TEnv*)inFile_p->Get("config");  

  const Int_t nMaxJetRs = 7;
  std::string jetRStr = inFileConfig_p->GetValue("JETRS", "");
  std::vector<int> jetRVect = strToVectI(jetRStr);
  if(jetRVect.size() <= 0 || jetRVect.size() > (unsigned int)nMaxJetRs){
    std::cout << "Jet R query returned \'" << jetRStr << "\', either 0 valid jet Rs or in excess of cap, " << nMaxJetRs << ". return 1" << std::endl;
    return 1;
  }
  
  //Add some parameters from infileconfig to config for writeout
  config_p->SetValue("CENT", inFileConfig_p->GetValue("CENT", ""));
  config_p->SetValue("JETRS", jetRStr.c_str());
  
  Double_t evtWeight_;
  
  const Int_t nMaxPart = 10000;
  Int_t nPart_;
  Float_t pt_[nMaxPart];
  Float_t eta_[nMaxPart];
  Float_t phi_[nMaxPart];
  Float_t m_[nMaxPart];
  Int_t pid_[nMaxPart];

  const Int_t nMaxJets = 100;
  Int_t nJt_[nMaxJetRs];
  Float_t jtpt_[nMaxJetRs][nMaxJets];
  Float_t jtphi_[nMaxJetRs][nMaxJets];
  Float_t jteta_[nMaxJetRs][nMaxJets];
  Float_t jtm_[nMaxJetRs][nMaxJets];

  jewelTree_p->SetBranchStatus("*", 0);  
  jewelTree_p->SetBranchStatus("evtWeight", 1);
  jewelTree_p->SetBranchStatus("nPart", 1);  
  jewelTree_p->SetBranchStatus("pt", 1);  
  jewelTree_p->SetBranchStatus("eta", 1);  
  jewelTree_p->SetBranchStatus("phi", 1);  
  jewelTree_p->SetBranchStatus("m", 1);  
  jewelTree_p->SetBranchStatus("pid", 1);  

  jewelTree_p->SetBranchAddress("evtWeight", &evtWeight_);
  jewelTree_p->SetBranchAddress("nPart", &nPart_);
  jewelTree_p->SetBranchAddress("pt", pt_);
  jewelTree_p->SetBranchAddress("eta", eta_);
  jewelTree_p->SetBranchAddress("phi", phi_);
  jewelTree_p->SetBranchAddress("m", m_);
  jewelTree_p->SetBranchAddress("pid", pid_);

  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    jewelTree_p->SetBranchStatus(("nJtR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtptR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtetaR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtphiR" + std::to_string(jetRVect[rI])).c_str(), 1);
    jewelTree_p->SetBranchStatus(("jtmR" + std::to_string(jetRVect[rI])).c_str(), 1);

    jewelTree_p->SetBranchAddress(("nJtR" + std::to_string(jetRVect[rI])).c_str(), &nJt_[rI]);
    jewelTree_p->SetBranchAddress(("jtptR" + std::to_string(jetRVect[rI])).c_str(), jtpt_[rI]);
    jewelTree_p->SetBranchAddress(("jtetaR" + std::to_string(jetRVect[rI])).c_str(), jteta_[rI]);
    jewelTree_p->SetBranchAddress(("jtphiR" + std::to_string(jetRVect[rI])).c_str(), jtphi_[rI]);
    jewelTree_p->SetBranchAddress(("jtmR" + std::to_string(jetRVect[rI])).c_str(), jtm_[rI]);
  }
  
  //Check the outrootfile directory exists or add directory path to outfilename
  if(outFileName.find("/") == std::string::npos){
    outFileName = outputStr + "/" + outFileName;
  }
  else{
    std::string dirName = outFileName.substr(0, outFileName.rfind("/"));
    check.doCheckMakeDir(dirName);
  }
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");

  //Construct the histograms, which can cover many R values and gamma pt bins
  ULong64_t nPhotonsPerPtBin[nMaxPtBins];
  TH1D* varHist_p[nMaxJetRs][nMaxPtBins];
  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
      nPhotonsPerPtBin[gI] = 0;

      std::string nameStr = varNameLower + "_R" + std::to_string(jetRVect[rI]) + "_GammaPt" + std::to_string(gI) + "_h";
      std::string titleStr = ";" + varNameLabel + ";#frac{1}{N_{#gamma}}";
      if(isMultijet) titleStr = titleStr + "#frac{dN_{JJ#gamma}}{d" + varNameLabel + "}";
      else titleStr = titleStr + "#frac{dN_{J#gamma}}{d" + varNameLabel + "}";
      
      varHist_p[rI][gI] = new TH1D(nameStr.c_str(), titleStr.c_str(), nVarBins, varBins);
      setSumW2(varHist_p[rI][gI]);
    }
  }
  
  const ULong64_t nEntries = jewelTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  std::cout << "Processing " << nEntries << " events..." << std::endl;
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry % nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;
    jewelTree_p->GetEntry(entry);

    std::vector<TLorentzVector> goodPhotons;
    std::vector<std::vector<TLorentzVector> > goodJets;
    
    //Process the event for every valid possible photon
    for(Int_t pI = 0; pI < nPart_; ++pI){
      if(pid_[pI] != 22) continue; //Only check photons
      //Gamma Pt Cuts
      if(pt_[pI] < gammaPtBinsLowReco) continue;
      if(pt_[pI] > gammaPtBinsHighReco) continue;
      //Gamma Eta Cuts
      if(!photonEtaIsGood(eta_[pI])) continue;

      //We need to impose isolation condition
      Float_t genEtSum4 = 0.0;
      for(Int_t pI2 = 0; pI2 < nPart_; ++pI2){
	//Skip the particle itself
	if(pI == pI2) continue;

	//Skip muons and all neutrinos
	if(TMath::Abs(pid_[pI2]) == 13) continue;
	if(TMath::Abs(pid_[pI2]) == 12) continue;
	if(TMath::Abs(pid_[pI2]) == 14) continue;
	if(TMath::Abs(pid_[pI2]) == 16) continue;

	//Only keep particles w/ dR < 0.4 (just hard-code it)
	Float_t dR = getDR(eta_[pI], phi_[pI], eta_[pI2], phi_[pI2]);
	if(dR > 0.4) continue;

	//Add to the sumEtive
	TLorentzVector temp;
	temp.SetPtEtaPhiM(pt_[pI2], eta_[pI2], phi_[pI2], m_[pI2]);

	genEtSum4 += temp.Et();
      }
      if(genEtSum4 > 5.0) continue;

      TLorentzVector tL;
      tL.SetPtEtaPhiM(pt_[pI], eta_[pI], phi_[pI], 0.0);
      goodPhotons.push_back(tL);
    }

    for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
      goodJets.push_back({});
      
      for(Int_t jI = 0; jI < nJt_[rI]; ++jI){
	//Jet Pt Cuts
	if(jtpt_[rI][jI] < jtPtBinsLowReco) continue;
	if(jtpt_[rI][jI] > jtPtBinsHighReco) continue;
	//Jet Eta Cuts
	if(jteta_[rI][jI] < jtEtaBinsLow) continue;
	if(jteta_[rI][jI] > jtEtaBinsHigh) continue;

	TLorentzVector tL;
	tL.SetPtEtaPhiM(jtpt_[rI][jI], jteta_[rI][jI], jtphi_[rI][jI], jtm_[rI][jI]);
	goodJets[rI].push_back(tL);
      }      
    }

    //Construct multi-jet observables
    for(unsigned int gI = 0; gI < goodPhotons.size(); ++gI){      
      Int_t gammaPos = ghostPos(gammaPtBinsVect, goodPhotons[gI].Pt());
      ++nPhotonsPerPtBin[gammaPos];
      
      for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
	for(unsigned int jI = 0; jI < goodJets[rI].size(); ++jI){
	  //Exclude jets w/ association w/ photon
	  Double_t dRJG = getDR(goodPhotons[gI].Eta(), goodPhotons[gI].Phi(), goodJets[rI][jI].Eta(), goodJets[rI][jI].Phi());
	  if(dRJG < gammaJtDRExclusionCut) continue;
	  
	  //DPhi between jet and gamma
	  Double_t dPhiJG = getDPHI(goodPhotons[gI].Phi(), goodJets[rI][jI].Phi());
	  if(TMath::Abs(dPhiJG) < gammaJtDPhiCut) continue;
	  
	  if(!isMultijet){
	    //Give it the same jet twice - in inclusive jets its not used
	    Float_t varVal = getVar(varNameLower, goodJets[rI][jI], goodJets[rI][jI], goodPhotons[gI]);
	    varHist_p[rI][gammaPos]->Fill(varVal, evtWeight_);
	  }
	  else{
	    //Loop over jets again
	    for(unsigned int jI2 = jI+1; jI2 < goodJets[rI].size(); ++jI2){
	      //Gotta impose more cuts; same cuts as on prev jet, but adding multijet cuts
	      Double_t dRJG2 = getDR(goodPhotons[gI].Eta(), goodPhotons[gI].Phi(), goodJets[rI][jI2].Eta(), goodJets[rI][jI2].Phi());
	      if(dRJG2 < gammaJtDRExclusionCut) continue;
	      
	      //DPhi between jet and gamma
	      Double_t dPhiJG2 = getDPHI(goodPhotons[gI].Phi(), goodJets[rI][jI2].Phi());
	      if(TMath::Abs(dPhiJG2) < gammaJtDPhiCut) continue;

	      //do drjj cut
	      Double_t dRJJ = getDR(goodJets[rI][jI].Eta(), goodJets[rI][jI].Phi(), goodJets[rI][jI2].Eta(), goodJets[rI][jI2].Phi());
	      if(dRJJ < mixJtDRExclusionCut) continue;
	      
	      //Now construct multijet and do multijet dphi cut
	      TLorentzVector multiJt = goodJets[rI][jI] + goodJets[rI][jI2];
	      Double_t dPhiJJG = getDPHI(goodPhotons[gI].Phi(), multiJt.Phi());
	      if(TMath::Abs(dPhiJJG) < gammaMultiJtDPhiCut) continue;

	      Float_t varVal = getVar(varNameLower, goodJets[rI][jI], goodJets[rI][jI2], goodPhotons[gI]);
	      varHist_p[rI][gammaPos]->Fill(varVal, evtWeight_);
	    }
	  }	  
	}	
      }      
    }    
  }
  std::cout << "Event processing complete." << std::endl;    
  
  inFile_p->Close();
  delete inFile_p;

  outFile_p->cd();


  for(unsigned int rI = 0; rI < jetRVect.size(); ++rI){
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){    
      //Skip the ones outside your range
      Double_t binCenter = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;
      if(binCenter < gammaPtBinsLowReco) continue;
      if(binCenter > gammaPtBinsHighReco) continue;
      
      //Scale by nPhotons
      binWidthAndScaleNorm(varHist_p[rI][gI], (Double_t)nPhotonsPerPtBin[gI]);

      varHist_p[rI][gI]->Write("", TObject::kOverwrite);
      delete varHist_p[rI][gI];
    }
  }
  
  config_p->Write("config", TObject::kOverwrite);
  
  outFile_p->Close();
  delete outFile_p;
  
  //Clean all new's
  delete config_p;  

  std::cout << "GDJHEPMCANA complete. return 0" << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHEPMCAna.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHEPMCAna(argv[1]);
  return retVal;
}
