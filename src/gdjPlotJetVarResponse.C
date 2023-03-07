//Author: Chris McGinn (2023.02.28)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs  

//c+cpp
#include <iostream>
#include <string>
#include <vector>

//ROOT
#include "TCanvas.h"
#include "TEnv.h"
#include "TFile.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TMath.h"

//Local
#include "include/binUtils.h"
#include "include/checkMakeDir.h"
#include "include/envUtil.h"
#include "include/globalDebugHandler.h"
#include "include/plotUtilities.h"
#include "include/stringUtil.h"


int gdjPlotJetVarResponse(const std::string inConfigFileName)
{
  //Check config exists
  checkMakeDir check;
  if(!check.checkFileExt(inConfigFileName, ".config")) return 1;

  //Create date-stamped pdfDir
  const std::string dateStr = getDateStr();
  const std::string pdfDirStr = "pdfDir/" + dateStr;
  check.doCheckMakeDir("pdfDir");
  check.doCheckMakeDir(pdfDirStr);

  //Global debug handler - do 'export DOGLOBALDEBUGROOT=1' to turn on debug output
  globalDebugHandler gDebug;
  const bool doGlobalDebug = gDebug.GetDoGlobalDebug();
  //Grab our input config, define necessary internal params
  TEnv* config_p = new TEnv(inConfigFileName.c_str());
  std::vector<std::string> necessaryParams = {
    "INUNFOLDFILENAME",
    "NLABEL",
    "LABELX.0",
    "LABELY.0",
    "MINZ",
    "MAXZ",
    "DOLOGZ"};
  

  //Check that the required parameters exist
  if(!checkEnvForParams(config_p, necessaryParams)) return 1;
  
  //Grab in file and check it exists
  const std::string inUnfoldFileName = config_p->GetValue("INUNFOLDFILENAME", "");
  if(!check.checkFileExt(inUnfoldFileName, ".root")) return 1;

  //Grab the nLabels, labelXs and labelYs
  const Int_t nLabel = config_p->GetValue("NLABEL", -1);
  if(nLabel <= 0){
    std::cout << "WARNING: NLABEL=" << nLabel << " is <= 0, at least 1 must be specified. return 1" << std::endl;
    return 1;    
  }

  std::vector<Double_t> labelXs, labelYs;
  for(Int_t lI = 0; lI < 100; ++lI){
    Double_t labelX = config_p->GetValue(("LABELX." + std::to_string(lI)).c_str(), -1.0);
    Double_t labelY = config_p->GetValue(("LABELY." + std::to_string(lI)).c_str(), -1.0);

    if(labelX < 0.0 || labelX > 1.0){
      //      std::cout << "WARNING: LABELX." << lI << "=" << labelX << " is outside permitted range 0.0-1.0, continuing" << std::endl;
      continue;
    }
    else if(labelY < 0.0 || labelY > 1.0){
      //      std::cout << "WARNING: LABELY." << lI << "=" << labelY << " is outside permitted range 0.0-1.0, continuing" << std::endl;
      continue;
    }

    labelXs.push_back(labelX);
    labelYs.push_back(labelY);
  }
  
  Bool_t doLabel = true;
  if(labelXs.size() != (unsigned int)nLabel){
    std::cout << "WARNING: NLABEL=" << nLabel << " does not match size of labelXs and labelYs, " << labelXs.size() << ", turning off dolabel" << std::endl;
    doLabel = false;
  }

  std::vector<std::string> globalLabels;
  if(doLabel){
    for(unsigned int i = 0; i < 10; ++i){
      std::string globalLabel = config_p->GetValue(("GLOBALLABEL." + std::to_string(i)).c_str(), "");
      if(globalLabel.size() == 0) continue;
      globalLabels.push_back(globalLabel);
    }
  }

  const Double_t minZ = config_p->GetValue("MINZ", 100.0);
  const Double_t maxZ = config_p->GetValue("MAXZ", -100.0);
  const Bool_t doLogZ = (Bool_t)config_p->GetValue("DOLOGZ", 0);
					    
  if(minZ > maxZ){
    std::cout << "gdjPlotJetVarResponse: MINZ \'" << minZ << "\' is greater than MAXZ \'" << maxZ << "\'. Please fix. return w/o plotting (1)" << std::endl;
    return 1;
  }

  if(doLogZ && minZ < TMath::Power(10, -100)){
    std::cout << "gdjPlotJetVarResponse: DOLOGZ is 'true' but minZ <= 0.0 (='" << minZ << "'). Please fix. return w/o plotting (1)" << std::endl;
    return 1;
  }
  
  if(doGlobalDebug) std::cout << "DOGLOBALDEBUG FILE " << __FILE__ << ", L" << __LINE__ << std::endl;
  
  //Grab in file and relevant config / params
  TFile* inUnfoldFile_p = new TFile(inUnfoldFileName.c_str(), "READ");
  TEnv* inUnfoldConfig_p = (TEnv*)inUnfoldFile_p->Get("config");

  std::vector<std::string> necessaryUnfoldParams = {
    "NGAMMAPTBINS",
    "GAMMAPTBINSLOW",
    "GAMMAPTBINSHIGH",
    "GAMMAPTBINSLOWRECO",
    "GAMMAPTBINSHIGHRECO",
    "GAMMAPTBINSDOLOG",
    "GAMMAPTBINSDOCUSTOM",
    "VARNAME",
    "ISPP"};

  //Check that the required parameters exist
  if(!checkEnvForParams(inUnfoldConfig_p, necessaryUnfoldParams)) return 1;

  if(doGlobalDebug) std::cout << "DOGLOBALDEBUG FILE " << __FILE__ << ", L" << __LINE__ << std::endl;
  
  const Int_t nMaxBins = 200;

  //GAMMAPT BINS handling
  const Int_t nGammaPtBins = inUnfoldConfig_p->GetValue("NGAMMAPTBINS", 0);
  const Float_t gammaPtBinsLow = inUnfoldConfig_p->GetValue("GAMMAPTBINSLOW", 80.0);
  const Float_t gammaPtBinsHigh = inUnfoldConfig_p->GetValue("GAMMAPTBINSHIGH", 280.0);
  const Float_t gammaPtBinsLowReco = inUnfoldConfig_p->GetValue("GAMMAPTBINSLOWRECO", 80.0);
  const Float_t gammaPtBinsHighReco = inUnfoldConfig_p->GetValue("GAMMAPTBINSHIGHRECO", 280.0);
  const Bool_t gammaPtBinsDoLog = (bool)inUnfoldConfig_p->GetValue("GAMMAPTBINSDOLOG", 1);
  const Bool_t gammaPtBinsDoCustom = (bool)inUnfoldConfig_p->GetValue("GAMMAPTBINSDOCUSTOM", 1);
  Double_t gammaPtBins[nMaxBins+1];
  customBinFiller("GAMMAPT", inUnfoldConfig_p, nMaxBins, gammaPtBins);

  const std::string varName = inUnfoldConfig_p->GetValue("VARNAME", "");

  const bool isPP = inUnfoldConfig_p->GetValue("ISPP", 0);
  std::vector<std::string> centBinsStr = {"PP"};
  std::vector<std::string> centBinsLabel = {"#bf{p+p}"};

  std::vector<std::string> pbpbParams = {"CENTBINS"};
  
  if(!isPP){
    if(!checkEnvForParams(inUnfoldConfig_p, pbpbParams)) return 1;
    std::vector<std::string> tempCentBinsStr = strToVect(inUnfoldConfig_p->GetValue("CENTBINS", ""));
    centBinsStr.clear();
    centBinsLabel.clear();
    for(unsigned int cI = 0; cI < tempCentBinsStr.size()-1; ++cI){
      centBinsStr.push_back("Cent" + tempCentBinsStr[cI] + "to" + tempCentBinsStr[cI+1]);
      centBinsLabel.push_back("#bf{Pb+Pb, " + tempCentBinsStr[cI] + "-" + tempCentBinsStr[cI+1] + "%}");
    }
  }

   
  const int titleFont = 42;
  const double titleSize = 0.03;
 
  //Process all gamma pt bins
  for(Int_t gI = 0; gI < nGammaPtBins; ++gI){
    Float_t gammaPtBinCenter = (gammaPtBins[gI] + gammaPtBins[gI+1])/2.0;
    if(gammaPtBinCenter < gammaPtBinsLowReco) continue;
    if(gammaPtBinCenter > gammaPtBinsHighReco) continue;

    for(unsigned int cI = 0; cI < centBinsStr.size(); ++cI){
      TCanvas* canv_p = new TCanvas("canv", "canv", 1100, 900);
      canv_p->SetRightMargin(canv_p->GetRightMargin()*1.2);
      
      std::string jetVarHistName = "rooResJt" + varName + "_JetVarMatrixForRecoPho_" + centBinsStr[cI] + "_GammaPt" + std::to_string(gI) + "_Nominal_h";
      TH2D* jetVarRes_p = (TH2D*)inUnfoldFile_p->Get(jetVarHistName.c_str());

      jetVarRes_p->Scale(1.0/jetVarRes_p->Integral());
      
      jetVarRes_p->SetMinimum(minZ);
      jetVarRes_p->SetMaximum(maxZ);
      jetVarRes_p->DrawCopy("COLZ");
      if(doLogZ) canv_p->SetLogz();
      
      if(doLabel){
	std::vector<std::string> currLabels = globalLabels;
	currLabels.push_back(prettyString(gammaPtBins[gI], 1, false) + " < p_{T,#gamma}^{Truth} < " + prettyString(gammaPtBins[gI+1], 1, false));
	currLabels.push_back(centBinsLabel[cI]);
	
	TLatex* label_p = new TLatex();
	initLabel(label_p, titleFont, titleSize, 0);
	
	//Quick calculation of how to do the labels
	Int_t currLabelsSize = currLabels.size();
	Int_t labelPerPos = currLabelsSize/nLabel;
	if(labelPerPos < 2) labelPerPos = 2;
	
	for(unsigned int cI = 0; cI < currLabels.size(); ++cI){
	  label_p->DrawLatex(labelXs[cI/labelPerPos], labelYs[cI/labelPerPos] - 0.045*(cI%labelPerPos), currLabels[cI].c_str());
	}
	delete label_p;
      }
      
      std::string saveName = pdfDirStr + "/jetVarResponse" + varName + "_GammaPt" + std::to_string(gI) + "_" + centBinsStr[cI]  + "_" + dateStr + ".png"; 
      quietSaveAs(canv_p, saveName);
      delete canv_p;
    }
  }

  if(doGlobalDebug) std::cout << "DOGLOBALDEBUG FILE " << __FILE__ << ", L" << __LINE__ << std::endl;
  
  inUnfoldFile_p->Close();
  
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/gdjPlotJetVarResponse.exe <inConfigFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }

  int retVal = 0;
  retVal += gdjPlotJetVarResponse(argv[1]);
  return retVal;
}
