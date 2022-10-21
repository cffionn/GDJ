
//c and cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPad.h"
#include "TStyle.h"
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
#include "include/HIJetPlotStyle.h"
#include "include/histDefUtility.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"
#include "include/varUtil.h"

bool isQuarkFromPDGID(int inPDGID)
{
  bool isQuark = false;
  if(TMath::Abs(inPDGID) == 1) isQuark = true;  
  else if(TMath::Abs(inPDGID) == 2) isQuark = true;  
  else if(TMath::Abs(inPDGID) == 3) isQuark = true;  
  else if(TMath::Abs(inPDGID) == 4) isQuark = true;  
  else if(TMath::Abs(inPDGID) == 5) isQuark = true;  
  else if(TMath::Abs(inPDGID) == 6) isQuark = true;    
  return isQuark;
}

bool isGluonFromPDGID(int inPDGID)
{
  bool isGluon = false;
  if(inPDGID == 21) isGluon = true;
  return isGluon;
}

int gdjHistToGenVarPlots(const std::string inConfigFileName)
{
  //Class for checking file/dir existence + creation
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;  

  //Debugging is set as an environment variable - we use a simple class to pick up
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Date string for tagging all output
  const std::string dateStr = getDateStr();
  const std::string outputStr = "output/" + dateStr;
  const std::string pdfStr = "pdfDir/" + dateStr;
  //output of current job will all go to date tagged subdir
  check.doCheckMakeDir("output");
  check.doCheckMakeDir(outputStr);

  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfStr);

  //Get config
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  //These params are needed to run - other params are optional
  std::vector<std::string> necessaryParams = {"INRESPONSEFILENAME",
					      "OUTFILENAME"};

  //Terminate if not all necessary params are found
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;

  std::string inResponseFileName = config_p->GetValue("INRESPONSEFILENAME", "");
  //check that our given input files exist
  if(!check.checkFileExt(inResponseFileName, ".root")) return 1;

  //if response file exists open, and grab internal config object
  TFile* inResponseFile_p = new TFile(inResponseFileName.c_str(), "READ");
  TEnv* inResponseFileConfig_p = (TEnv*)inResponseFile_p->Get("config");
  
  std::vector<std::string> necessaryInFileParams = {"ISMC",
						    "NGAMMAPTBINS",
                                                    "GAMMAPTBINSLOW",
                                                    "GAMMAPTBINSHIGH",
                                                    "GAMMAPTBINSDOLOG",
                                                    "GAMMAPTBINSDOCUSTOM",
						    "NJTPTBINS",
                                                    "JTPTBINSLOW",
                                                    "JTPTBINSHIGH",
                                                    "JTPTBINSDOLOG",
                                                    "JTPTBINSDOCUSTOM",
						    "JTETABINSLOW",
						    "JTETABINSHIGH",
						    "GAMMAJTDPHI",
						    "GAMMAMULTIJTDPHI",
						    "MIXJETEXCLUSIONDR"};

  //Return failed if we are missing required parameters
  if(!checkEnvForParams(inResponseFileConfig_p, necessaryInFileParams)) return 1;

  const Int_t nMaxBins = 500;
  
  //Grab gamma and jt pt bins since they will define things for all plots
  //First grab the full set of gamma pt bins
  Int_t nGammaPtBins = inResponseFileConfig_p->GetValue("NGAMMAPTBINS", 0);
  Float_t gammaPtBinsLow = inResponseFileConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  Float_t gammaPtBinsHigh = inResponseFileConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  Bool_t gammaPtBinsDoLog = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  Bool_t gammaPtBinsDoCustom = (bool)inResponseFileConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxBins+1];

  if(gammaPtBinsDoLog) getLogBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);
  else if(gammaPtBinsDoCustom){
    std::string gammaPtBinsStr = inResponseFileConfig_p->GetValue("GAMMAPTBINSCUSTOM", "");
    if(gammaPtBinsStr.size() == 0){
      std::cout << "gammaPtBinsDoCustom is true but gammaPtBinsCustom is empty. return 1" << std::endl;
      return 1;
    }

    std::vector<float> gammaPtBinsTemp = strToVectF(gammaPtBinsStr);
    Int_t nGammaPtBinsTemp = ((Int_t)gammaPtBinsTemp.size()) - 1;
    if(nGammaPtBinsTemp != nGammaPtBins){
      std::cout << "Number of gamma pt custom bins doesnt match specified number; return 1" << std::endl;
      return 1;
    }

    for(unsigned int gI = 0; gI < gammaPtBinsTemp.size(); ++gI){
      gammaPtBins[gI] = gammaPtBinsTemp[gI];
    }
  }
  else getLinBins(gammaPtBinsLow, gammaPtBinsHigh, nGammaPtBins, gammaPtBins);

  //For jet bins we just need the max and min val for cutting
  const Float_t jtPtBinsLow = inResponseFileConfig_p->GetValue("JTPTBINSLOW", -1.0);
  const Float_t jtPtBinsHigh = inResponseFileConfig_p->GetValue("JTPTBINSHIGH", -1.0);
  const Float_t jtEtaBinsLow = inResponseFileConfig_p->GetValue("JTETABINSLOW", -1.0);
  const Float_t jtEtaBinsHigh = inResponseFileConfig_p->GetValue("JTETABINSHIGH", -1.0);

  const Float_t gammaJtDPhiCut = mathStringToNum(inResponseFileConfig_p->GetValue("GAMMAJTDPHI", ""));
  const Float_t gammaMultiJtDPhiCut = mathStringToNum(inResponseFileConfig_p->GetValue("GAMMAMULTIJTDPHI", ""));

  const Float_t mixJetExclusionDR = inResponseFileConfig_p->GetValue("MIXJETEXCLUSIONDR", -1.0);  
  
  std::vector<std::string> varNames = {"ajj", "xjj", "drjj"};
  std::vector<bool> goodVar;

  std::vector<std::string> qgVect = {"ALL", "QQ", "QG", "GQ", "GG"};
  std::vector<std::string> qgLabel = {"All Combinations", "q+q", "q+g", "g+q", "g+g"};
  
  //Macro index is gammaPt, micro index is var
  std::vector<std::vector<std::vector<TH1D*> > > histsPerGammaPtPerVarFlav;
  std::vector<Double_t> phoCountPerGammaPt;
  
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    histsPerGammaPtPerVarFlav.push_back({});
    phoCountPerGammaPt.push_back(0.0);

    for(unsigned int qgI = 0; qgI < qgVect.size(); ++qgI){
      histsPerGammaPtPerVarFlav[gI].push_back({});
    }
  }    

  std::string outFileName = config_p->GetValue("OUTFILENAME", "");
  outFileName = outputStr + "/" + outFileName;
  TFile* outFile_p = new TFile(outFileName.c_str(), "RECREATE");
  
  for(auto const& varName : varNames){
    std::cout << "VARNAME: " << varName << std::endl;
    std::string varNameUpper = returnAllCapsString(varName);

    if(isStrSame(varNameUpper, "DRJJ")){
      varNameUpper = "DR";
    }
    else if(isStrSame(varNameUpper, "AJJ")){
      varNameUpper = "AJ";
    }
    
    Int_t nVarBins = inResponseFileConfig_p->GetValue(("N" + varNameUpper + "BINS").c_str(), -1);
    //Check that the number of bins is sensible and within limits
    if(nVarBins <= 0 || nVarBins > nMaxBins){
      std::cout << "n" << varNameUpper << "Bin=" << nVarBins << " is not valid - please choose a value between 0-" << nMaxBins << ". continue on processing var \'" << varName << "\'." << std::endl;
      goodVar.push_back(false);
    }
    else goodVar.push_back(true);

    Double_t varBinsLow = inResponseFileConfig_p->GetValue((varNameUpper + "BINSLOW").c_str(), -1.0);
    Double_t varBinsHigh = inResponseFileConfig_p->GetValue((varNameUpper + "BINSHIGH").c_str(), -1.0);
    Bool_t varBinsDoLog = (Bool_t)inResponseFileConfig_p->GetValue((varNameUpper + "BINSDOLOG").c_str(), 0);
    Bool_t varBinsDoCustom = (Bool_t)inResponseFileConfig_p->GetValue((varNameUpper + "BINSDOCUSTOM").c_str(), 0);
    Double_t varBins[nMaxBins+1];

    if(varBinsDoLog) getLogBins(varBinsLow, varBinsHigh, nVarBins, varBins);
    else if(varBinsDoCustom){
      std::string binsCustomStr = inResponseFileConfig_p->GetValue((varNameUpper + "BINSCUSTOM").c_str(), "");
      std::vector<float> binsCustomVect = strToVectF(binsCustomStr);
      for(unsigned int i = 0; i < binsCustomVect.size(); ++i){
	varBins[i] = binsCustomVect[i];
      }
    }
    else getLinBins(varBinsLow, varBinsHigh, nVarBins, varBins);
    
    std::cout << " ";
    for(Int_t bI = 0; bI < nVarBins; ++bI){
      std::cout << varBins[bI] << ", ";
    }
    std::cout << varBins[nVarBins] << "." << std::endl;

    std::string varLabel = varNameToLabel(varName);
    std::string titleName = ";" + varLabel + ";#frac{1}{N_{#gamma}} #frac{d" + varNameToNLabel(varName) + "}{d" + varLabel + "}";
    for(Int_t gI = 0; gI < nGammaPtBins; ++gI){            
      for(unsigned int qgI = 0; qgI < qgVect.size(); ++qgI){	
	std::string histName = varName + "Gen_GammaPt" + std::to_string(gI) + "_" + qgVect[qgI] + "_h";
	
	histsPerGammaPtPerVarFlav[gI][qgI].push_back(new TH1D(histName.c_str(), titleName.c_str(), nVarBins, varBins));
	histsPerGammaPtPerVarFlav[gI][qgI][histsPerGammaPtPerVarFlav[gI][qgI].size()-1]->Sumw2();
      }
    }
  }    

  bool oneGoodVar = false;
  for(unsigned int gI = 0; gI < goodVar.size(); ++gI){
    oneGoodVar = oneGoodVar || goodVar[gI];
  }

  if(!oneGoodVar){
    std::cout << "No variable is found. return 1" << std::endl;
    return 1;
  }


  //Now start processing by defining the set of variables we need from the ttree
  Float_t truthGammaPt_;
  Float_t truthGammaPhi_;
  Float_t truthGammaEta_;

  const Int_t nMaxJets = 100;
  Int_t nRecoJt_;
  Float_t truthJtPt_[nMaxJets];
  Float_t truthJtPhi_[nMaxJets];
  Float_t truthJtEta_[nMaxJets];
  Int_t truthJtFlavor_[nMaxJets];

  Int_t nTruthJtUnmatched_;
  Float_t truthJtUnmatchedPt_[nMaxJets];
  Float_t truthJtUnmatchedPhi_[nMaxJets];
  Float_t truthJtUnmatchedEta_[nMaxJets];
  Int_t truthJtUnmatchedFlavor_[nMaxJets];

  Double_t unfoldWeight_; 

  //Grab the unfolding ttree for hist filling
  TTree* unfoldTree_p = (TTree*)inResponseFile_p->Get("unfoldTree_p");
  unfoldTree_p->SetBranchStatus("*", 0);

  unfoldTree_p->SetBranchStatus("truthGammaPt", 1);
  unfoldTree_p->SetBranchStatus("truthGammaEta", 1);
  unfoldTree_p->SetBranchStatus("truthGammaPhi", 1);
  
  unfoldTree_p->SetBranchStatus("nRecoJt", 1);
  unfoldTree_p->SetBranchStatus("truthJtPt", 1);
  unfoldTree_p->SetBranchStatus("truthJtEta", 1);
  unfoldTree_p->SetBranchStatus("truthJtPhi", 1);
  unfoldTree_p->SetBranchStatus("truthJtFlavor", 1);

  unfoldTree_p->SetBranchStatus("nTruthJtUnmatched", 1);
  unfoldTree_p->SetBranchStatus("truthJtUnmatchedPt", 1);
  unfoldTree_p->SetBranchStatus("truthJtUnmatchedEta", 1);
  unfoldTree_p->SetBranchStatus("truthJtUnmatchedPhi", 1);
  unfoldTree_p->SetBranchStatus("truthJtUnmatchedFlavor", 1);

  unfoldTree_p->SetBranchStatus("unfoldWeight", 1);
  
  unfoldTree_p->SetBranchAddress("truthGammaPt", &truthGammaPt_);
  unfoldTree_p->SetBranchAddress("truthGammaEta", &truthGammaEta_);
  unfoldTree_p->SetBranchAddress("truthGammaPhi", &truthGammaPhi_);
  
  unfoldTree_p->SetBranchAddress("nRecoJt", &nRecoJt_);
  unfoldTree_p->SetBranchAddress("truthJtPt", truthJtPt_);
  unfoldTree_p->SetBranchAddress("truthJtPhi", truthJtPhi_);
  unfoldTree_p->SetBranchAddress("truthJtEta", truthJtEta_);
  unfoldTree_p->SetBranchAddress("truthJtFlavor", truthJtFlavor_);

  unfoldTree_p->SetBranchAddress("nTruthJtUnmatched", &nTruthJtUnmatched_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPt", truthJtUnmatchedPt_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedPhi", truthJtUnmatchedPhi_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedEta", truthJtUnmatchedEta_);
  unfoldTree_p->SetBranchAddress("truthJtUnmatchedFlavor", truthJtUnmatchedFlavor_);

  unfoldTree_p->SetBranchAddress("unfoldWeight", &unfoldWeight_);
  
  const ULong64_t nEntries = unfoldTree_p->GetEntries();
  const ULong64_t nDiv = TMath::Max((ULong64_t)1, nEntries/20);

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  for(ULong64_t entry = 0; entry < nEntries; ++entry){
    if(entry%nDiv == 0) std::cout << " Entry " << entry << "/" << nEntries << "..." << std::endl;    
    unfoldTree_p->GetEntry(entry);

    //Grab photon position
    Int_t gammaPtPos = ghostPos(nGammaPtBins, gammaPtBins, truthGammaPt_);
    if(gammaPtPos < 0) continue;
    if(!photonEtaIsGood(truthGammaEta_)) continue;

    phoCountPerGammaPt[gammaPtPos] += unfoldWeight_;
    
    TLorentzVector phoTL;
    phoTL.SetPtEtaPhiM(truthGammaPt_, truthGammaEta_, truthGammaPhi_, 0.0);
    
    //Get all truth jets into a single vector
    std::vector<TLorentzVector> truthJets;
    std::vector<Int_t> truthFlavors;
    for(Int_t tI = 0; tI < nRecoJt_; ++tI){
      if(truthJtPt_[tI] < jtPtBinsLow) continue;
      if(truthJtPt_[tI] >= jtPtBinsHigh) continue;
      if(truthJtEta_[tI] <= jtEtaBinsLow) continue;
      if(truthJtEta_[tI] >= jtEtaBinsHigh) continue;
      
      Float_t jtDPhi = TMath::Abs(getDPHI(truthJtPhi_[tI], truthGammaPhi_));
      if(jtDPhi < gammaJtDPhiCut) continue;
      
      TLorentzVector tempTL;
      tempTL.SetPtEtaPhiM(truthJtPt_[tI], truthJtEta_[tI], truthJtPhi_[tI], 0.0);
      truthJets.push_back(tempTL);
      truthFlavors.push_back(truthJtFlavor_[tI]);
    }

    for(Int_t tI = 0; tI < nTruthJtUnmatched_; ++tI){
      if(truthJtUnmatchedPt_[tI] < jtPtBinsLow) continue;
      if(truthJtUnmatchedPt_[tI] >= jtPtBinsHigh) continue;
      if(truthJtUnmatchedEta_[tI] <= jtEtaBinsLow) continue;
      if(truthJtUnmatchedEta_[tI] >= jtEtaBinsHigh) continue;

      Float_t jtDPhi = TMath::Abs(getDPHI(truthJtUnmatchedPhi_[tI], truthGammaPhi_));
      if(jtDPhi < gammaJtDPhiCut) continue;

      TLorentzVector tempTL;
      tempTL.SetPtEtaPhiM(truthJtUnmatchedPt_[tI], truthJtUnmatchedEta_[tI], truthJtUnmatchedPhi_[tI], 0.0);
      truthJets.push_back(tempTL);      
      truthFlavors.push_back(truthJtUnmatchedFlavor_[tI]);
    }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
    
    for(unsigned int tI = 0; tI < truthJets.size(); ++tI){ 
      for(unsigned int tI2 = tI+1; tI2 < truthJets.size(); ++tI2){
	Float_t dR = getDR(truthJets[tI].Eta(), truthJets[tI].Phi(), truthJets[tI2].Eta(), truthJets[tI2].Phi());

	if(dR < mixJetExclusionDR) continue;
	
	TLorentzVector combinedJet = truthJets[tI] + truthJets[tI2];
	Float_t multiJtDPhi = TMath::Abs(getDPHI(combinedJet.Phi(), truthGammaPhi_));

	if(multiJtDPhi < gammaMultiJtDPhiCut) continue;
      
	std::vector<int> qgPos = {0};
	Int_t leadingFlavor = truthFlavors[tI];
	Int_t subleadingFlavor = truthFlavors[tI2];	
	if(truthJets[tI].Pt() < truthJets[tI2].Pt()){
	  leadingFlavor = truthFlavors[tI2];	  
	  subleadingFlavor = truthFlavors[tI];	  
	}	     	

	if(isQuarkFromPDGID(leadingFlavor) && isQuarkFromPDGID(subleadingFlavor)) qgPos.push_back(1);
	else if(isQuarkFromPDGID(leadingFlavor) && isGluonFromPDGID(subleadingFlavor)) qgPos.push_back(2);
	else if(isGluonFromPDGID(leadingFlavor) && isQuarkFromPDGID(subleadingFlavor)) qgPos.push_back(3);
	else if(isGluonFromPDGID(leadingFlavor) && isGluonFromPDGID(subleadingFlavor)) qgPos.push_back(4);
	
	
	for(unsigned int vI = 0; vI < varNames.size(); ++vI){
	  Float_t varVal = getVar(varNames[vI], truthJets[tI], truthJets[tI2], phoTL);

	  if(gammaPtPos == 3 && false){
	    if(isStrSame(varNames[vI], "drjj")){
	      if(varVal > 2.8){
		std::cout << "HIGH VAR VAL: " << varVal << std::endl;
		std::cout << " Pho pt eta phi: " << truthGammaPt_ << ", " << truthGammaEta_ << ", " << truthGammaPhi_ << std::endl;
		std::cout << " Jet 1 pt Eta Phi: " << truthJets[tI].Pt() << ", " << truthJets[tI].Eta() << ", " << truthJets[tI].Phi() << std::endl;
		std::cout << " Jet 2 pt Eta Phi: " << truthJets[tI2].Pt() << ", " << truthJets[tI2].Eta() << ", " << truthJets[tI2].Phi() << std::endl;
	      }
	    }
	  }

	  for(unsigned int qgI = 0; qgI < qgPos.size(); ++qgI){	      
	    histsPerGammaPtPerVarFlav[gammaPtPos][qgPos[qgI]][vI]->Fill(varVal, unfoldWeight_);
	  }	  
	}   
      }     
    }     
  }

  if(doGlobalDebug) std::cout << "FILE, LINE: " << __FILE__ << ", " << __LINE__ << std::endl;
  
  //Save the hists
  outFile_p->cd();

  for(int gI = 0; gI < nGammaPtBins; ++gI){
    for(unsigned int qgI = 0; qgI < qgVect.size(); ++qgI){	
      for(unsigned int vI = 0; vI < varNames.size(); ++vI){
	binWidthAndScaleNorm(histsPerGammaPtPerVarFlav[gI][qgI][vI], phoCountPerGammaPt[gI]);
	histsPerGammaPtPerVarFlav[gI][qgI][vI]->Write("", TObject::kOverwrite);
      }
    }    
  }

  //Define some variables for plotting
  const Double_t leftMargin = 0.15;
  const Double_t bottomAbsFracOfLeft = 0.7;
  const Double_t bottomMargin = bottomAbsFracOfLeft*leftMargin;
  const Double_t rightMargin = 0.03;
  const Double_t topMargin = rightMargin;

  const int titleFont = 42;
  const double titleSize = 0.035;

  //Make some plots
  for(int gI = 0; gI < nGammaPtBins; ++gI){
    for(unsigned int vI = 0; vI < varNames.size(); ++vI){

      TCanvas* canv_p = new TCanvas("canv_p", "", 900, 900);
      setMargins(canv_p, leftMargin, topMargin, rightMargin, leftMargin);

      HIJet::Style::EquipHistogram(histsPerGammaPtPerVarFlav[gI][0][vI], 0);
      histsPerGammaPtPerVarFlav[gI][0][vI]->SetMarkerSize(0.001);
      std::cout << "OFFSET: " << histsPerGammaPtPerVarFlav[gI][0][vI]->GetYaxis()->GetTitleOffset() << std::endl;
      histsPerGammaPtPerVarFlav[gI][0][vI]->GetYaxis()->SetTitleOffset(1.6);

      histsPerGammaPtPerVarFlav[gI][0][vI]->DrawCopy("HIST E1 P");

      gPad->SetTicks();
      
      TLegend* leg_p = new TLegend(0.6, 0.6, 0.9, 0.9);
      leg_p->SetBorderSize(0);
      leg_p->SetFillColor(0);
      leg_p->SetTextFont(titleFont);
      leg_p->SetTextSize(titleSize);
      
      for(unsigned int qgI = 1; qgI < qgVect.size(); ++qgI){	
	HIJet::Style::EquipHistogram(histsPerGammaPtPerVarFlav[gI][qgI][vI], qgI);
	histsPerGammaPtPerVarFlav[gI][qgI][vI]->SetFillColor(histsPerGammaPtPerVarFlav[gI][qgI][vI]->GetMarkerColor());

	TH1D* histTemp_p = (TH1D*)histsPerGammaPtPerVarFlav[gI][qgI][vI]->Clone("histTemp_h");

	for(unsigned int qgI2 = qgVect.size()-1; qgI2 > qgI; --qgI2){
	  histTemp_p->Add(histsPerGammaPtPerVarFlav[gI][qgI2][vI]);
	}
	
	HIJet::Style::EquipHistogram(histTemp_p, qgI);
	histTemp_p->SetFillColor(histTemp_p->GetMarkerColor());
	histTemp_p->DrawCopy("SAME HIST");
      }

      histsPerGammaPtPerVarFlav[gI][0][vI]->DrawCopy("HIST E1 P SAME");

      canv_p->RedrawAxis();
      
      for(unsigned int qgI = 0; qgI < qgVect.size(); ++qgI){
	leg_p->AddEntry(histsPerGammaPtPerVarFlav[gI][qgI][vI], qgLabel[qgI].c_str(), "P L");
      }

      leg_p->Draw("SAME");
      
      gStyle->SetOptStat(0);
      
      std::string saveName = pdfStr + "/truth" + returnAllCapsString(varNames[vI]) + "_GammaPt" + std::to_string(gI) + "_" + dateStr + ".png";
      quietSaveAs(canv_p, saveName);
      delete canv_p;
    }
  }
  
  //Cleanup the hists
  for(int gI = 0; gI < nGammaPtBins; ++gI){
    for(unsigned int qgI = 0; qgI < qgVect.size(); ++qgI){	
      for(unsigned int vI = 0; vI < varNames.size(); ++vI){
	delete histsPerGammaPtPerVarFlav[gI][qgI][vI];
      }
    }    
  }

  outFile_p->Close();
  delete outFile_p;
  
  std::cout << "GDJHISTTOUNFOLD COMPLETE. return 0." << std::endl;
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjHistToGenVarPlots.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjHistToGenVarPlots(argv[1]);
  return retVal;
}

