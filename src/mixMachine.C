//cpp dependencies
#include <iostream>

//ROOT dependencies
#include "THashList.h"

//Local dependencies
#include "include/envUtil.h"
#include "include/mixMachine.h"
#include "include/plotUtilities.h"

mixMachine::mixMachine(std::string inMixMachineName, mixMachine::mixMode inMixMode, TEnv* inParams_p)
{
  if(!Init(inMixMachineName, inMixMode, inParams_p)) std::cout << "MIXMACHINE: Initialization failure. return" << std::endl;
  return;
}

bool mixMachine::Init(std::string inMixMachineName, mixMachine::mixMode inMixMode, TEnv* inParams_p)
{
  m_mixMachineName = inMixMachineName;
  m_mixMode = inMixMode;

  std::vector<std::string> mixNames;
  if(m_mixMode == NONE) mixNames = m_noneMixNames;
  else if(m_mixMode == INCLUSIVE) mixNames = m_inclusiveMixNames;
  else if(m_mixMode == MULTI) mixNames = m_multiMixNames;
  
  THashList* paramList_p = (THashList*)inParams_p->GetTable();
  for(Int_t entry = 0; entry < paramList_p->GetEntries(); ++entry){
    std::string name = paramList_p->At(entry)->GetName();

    m_env.SetValue(name.c_str(), inParams_p->GetValue(name.c_str(), ""));
  }
  

  std::vector<std::string> reqParamsGlobal = {"IS2DUNFOLD",
					      "ISMC",
					      "NBINSX",
					      "BINSX",
					      "TITLEX"};

  std::vector<std::string> reqParams2DUnfold = {"NBINSY",
						"BINSY",
						"TITLEY"};
  
  if(!checkEnvForParams(&m_env, reqParamsGlobal)){
    Clean();
    return false;
  }

  m_is2DUnfold = m_env.GetValue("IS2DUNFOLD", 0);
  m_isMC = m_env.GetValue("ISMC", 0);
  m_nBinsX = m_env.GetValue("NBINSX", -1);
  m_binsXVect.clear();
  std::vector<Float_t> tempBins = strToVectF(m_env.GetValue("BINSX", ""));
  for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
    m_binsX[tI] = tempBins[tI];
    m_binsXVect.push_back(tempBins[tI]);
  }
  m_titleX = m_env.GetValue("TITLEX", "");  

  if(m_is2DUnfold){
    if(!checkEnvForParams(&m_env, reqParams2DUnfold)){
      Clean();
      return false;
    }

    m_nBinsY = m_env.GetValue("NBINSY", -1);
    m_binsYVect.clear();
    std::vector<Float_t> tempBins = strToVectF(m_env.GetValue("BINSY", ""));
    for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
      m_binsY[tI] = tempBins[tI];
      m_binsYVect.push_back(tempBins[tI]);
    }
    m_titleY = m_env.GetValue("TITLEY", "");  
    
    for(unsigned int mI = 0; mI < mixNames.size(); ++mI){
      if(!m_isMC && mixNames[mI].find("TRUTH") != std::string::npos) continue;

      m_hists2D.push_back(nullptr);

      std::string name = m_mixMachineName + "_MIXMODE" + std::to_string((int)m_mixMode) + "_" + mixNames[mI] + "_h";
      std::string title = ";" + m_titleX + ";" + m_titleY;      
      m_hists2D[mI] = new TH2F(name.c_str(), title.c_str(), m_nBinsX, m_binsX, m_nBinsY, m_binsY);
      m_hists2D[mI]->Sumw2();
    }
  }
  else{
    for(unsigned int mI = 0; mI < mixNames.size(); ++mI){
      if(!m_isMC && mixNames[mI].find("TRUTH") != std::string::npos) continue;

      m_hists1D.push_back(nullptr);

      std::string name = m_mixMachineName + "_MIXMODE" + std::to_string((int)m_mixMode) + "_" + mixNames[mI] + "_h";
      std::string title = ";" + m_titleX + ";Counts (Weighted)";      
      m_hists1D[mI] = new TH1F(name.c_str(), title.c_str(), m_nBinsX, m_binsX);
      m_hists1D[mI]->Sumw2();
    }    
  }

  m_isInit = true;
  return m_isInit;
}

bool mixMachine::FillXY(Float_t fillX, Float_t fillY, Float_t fillWeight, std::string mixName)
{
  if(!m_isInit){
    std::cout << "mixMachine::FillXY(): mixMachine is not initialized. return false" << std::endl;
  }
  
  int histPos = -1;

  std::string mixModeStr;
  std::vector<std::string> mixNames;

  if(m_mixMode == NONE){
    mixModeStr = "None";
    mixNames = m_noneMixNames;
  }
  else if(m_mixMode == INCLUSIVE){
    mixModeStr = "Inclusive";
    mixNames = m_inclusiveMixNames;
  }
  else if(m_mixMode == MULTI){
    mixModeStr = "Multi";
    mixNames = m_multiMixNames;
  }

  histPos = vectContainsStrPos(mixName, &mixNames);
  if(histPos < 0){
    std::cout << "mixMachine::FillXY(): Given mixName \'" << mixName << "\' not found for mixMode \'" << mixModeStr << "\'. Options are:" << std::endl;
    for(unsigned int i = 0; i < mixNames.size(); ++i){
      std::cout << " " << mixNames[i] << std::endl;
    }
    std::cout << "Returning false." << std::endl;
    return false;
  }  

  if(m_is2DUnfold){
    m_hists2D[histPos]->Fill(fillX, fillY, fillWeight);
  }
  else m_hists1D[histPos]->Fill(fillX, fillWeight);

  return true;
}

bool mixMachine::FillXYRaw(Float_t fillX, Float_t fillY, Float_t fillWeight)
{
  return FillXY(fillX, fillY, fillWeight, "RAW");
}

bool mixMachine::FillXYMix(Float_t fillX, Float_t fillY, Float_t fillWeight)
{
  return FillXY(fillX, fillY, fillWeight, "MIX");
}

bool mixMachine::FillXYMixCorrection(Float_t fillX, Float_t fillY, Float_t fillWeight)
{
  return FillXY(fillX, fillY, fillWeight, "MIXCORRECTION");
}

bool mixMachine::FillXYTruth(Float_t fillX, Float_t fillY, Float_t fillWeight)
{
  return FillXY(fillX, fillY, fillWeight, "TRUTH");
}

bool mixMachine::FillXYTruthMatchedReco(Float_t fillX, Float_t fillY, Float_t fillWeight)
{
  return FillXY(fillX, fillY, fillWeight, "TRUTHMATCHEDRECO");
}

bool mixMachine::FillX(Float_t fillX, Float_t fillWeight, std::string mixName)
{
  if(m_is2DUnfold){
    std::cout << "mixMachine::FillX(): Called despited m_is2DUnfold being true. return false" << std::endl;
    return false;
  }

  return FillXY(fillX, -1, fillWeight, mixName);
}

bool mixMachine::FillXRaw(Float_t fillX, Float_t fillWeight)
{
  return FillX(fillX, fillWeight, "RAW");
}

bool mixMachine::FillXMix(Float_t fillX, Float_t fillWeight)
{
  return FillX(fillX, fillWeight, "MIX");
}

bool mixMachine::FillXMixCorrection(Float_t fillX, Float_t fillWeight)
{
  return FillX(fillX, fillWeight, "MIXCORRECTION");
}

bool mixMachine::FillXTruth(Float_t fillX, Float_t fillWeight)
{
  return FillX(fillX, fillWeight, "TRUTH");
}

bool mixMachine::FillXTruthMatchedReco(Float_t fillX, Float_t fillWeight)
{
  return FillX(fillX, fillWeight, "TRUTHMATCHEDRECO");
}


bool mixMachine::CheckMachinesMatchMode(mixMachine* machineToCheck)
{
  mixMode checkMixMode = machineToCheck->GetMixMode();
  std::string machineToCheckName = machineToCheck->GetMixMachineName();

  if(checkMixMode != m_mixMode){
    std::cout << "MIXMACHINE::CHECKMACHINEMATCHMODE() ERROR: mixModes differ." << std::endl;
    std::cout << " Main Machine: " << m_mixMachineName << ", Mode: " << m_mixMode << std::endl;
    std::cout << " Machine-to-check: " << machineToCheckName << ", Mode: " << checkMixMode << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }
  return true;  
}


bool mixMachine::CheckMachinesMatchMC(mixMachine* machineToCheck)
{
  bool checkIsMC = machineToCheck->GetIsMC();
  std::string machineToCheckName = machineToCheck->GetMixMachineName();

  if(checkIsMC != m_isMC){
    std::cout << "MIXMACHINE::CHECKMACHINEMATCHMC() ERROR: isMC differs." << std::endl;
    std::cout << " Main Machine: " << m_mixMachineName << ", isMC: " << m_isMC << std::endl;
    std::cout << " Machine-to-check: " << machineToCheckName << ", isMC: " << checkIsMC << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  return true;
}

bool mixMachine::CheckMachinesMatch2D(mixMachine* machineToCheck)
{
  bool checkIs2D = machineToCheck->GetIs2DUnfold();
  std::string machineToCheckName = machineToCheck->GetMixMachineName();

  if(checkIs2D != m_is2DUnfold){
    std::cout << "MIXMACHINE::CHECKMACHINEMATCH2D() ERROR: is2DUnfold differs." << std::endl;
    std::cout << " Main Machine: " << m_mixMachineName << ", is2DUnfold: " << m_is2DUnfold << std::endl;
    std::cout << " Machine-to-check: " << machineToCheckName << ", isMC: " << checkIs2D << std::endl;
    std::cout << "return false" << std::endl;
    return false;
  }

  return true;
}

bool mixMachine::CheckMachinesMatchBins(mixMachine* machineToCheck, double precision)
{
  if(!this->CheckMachinesMatch2D(machineToCheck)){
    std::cout << "ABOVE ERROR COMES FROM INVOCATION OF MIXMACHINE::CHECKMACHINESMATCHBINS()" << std::endl;
    return false;
  }

  std::string machineToCheckName = machineToCheck->GetMixMachineName();

  if(m_is2DUnfold){
    Int_t nBinsXToCheck = machineToCheck->GetNBinsX();
    Int_t nBinsYToCheck = machineToCheck->GetNBinsY();

    if(nBinsXToCheck != m_nBinsX){
      std::cout << "MIXMACHINE::CHECKMACHINEMATCHBINS() ERROR: nBinsX differs." << std::endl;
      std::cout << " Main Machine: " << m_mixMachineName << ", nBinsX: " << m_nBinsX << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << ", nBinsX: "  << nBinsXToCheck << std::endl;
      std::cout << "return false" << std::endl;
      return false;
    }

    if(nBinsYToCheck != m_nBinsY){
      std::cout << "MIXMACHINE::CHECKMACHINEMATCHBINS() ERROR: nBinsY differs." << std::endl;
      std::cout << " Main Machine: " << m_mixMachineName << ", nBinsY: " << m_nBinsY << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << ", nBinsY: "  << nBinsYToCheck << std::endl;
      std::cout << "return false" << std::endl;
      return false;
    }

    std::vector<double> checkBinsX = machineToCheck->GetBinsX();
    std::vector<double> checkBinsY = machineToCheck->GetBinsY();

    bool allXBinsGood = true;
    for(unsigned int bIX = 0; bIX < checkBinsX.size(); ++bIX){
      if(TMath::Abs(checkBinsX[bIX] - m_binsXVect[bIX]) > precision){
	allXBinsGood = false;
	break;
      }
    }
    
    if(!allXBinsGood){
      std::cout << "MIXMACHINE::CheckMachinesMatchBins(): X-Binning mismatch." << std::endl;
      std::cout << " Precision level: " << precision << std::endl;
      std::cout << " Main machine: " << m_mixMachineName << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << std::endl;

      for(unsigned int bIX = 0; bIX < checkBinsX.size(); ++bIX){
	std::cout << " " << bIX << ": " << m_binsXVect[bIX] << ", " << checkBinsX[bIX] << std::endl;
      }
      
      return false;
    }

    bool allYBinsGood = true;
    for(unsigned int bIY = 0; bIY < checkBinsY.size(); ++bIY){
      if(TMath::Abs(checkBinsY[bIY] - m_binsYVect[bIY]) > precision){
	allYBinsGood = false;
	break;
      }
    }
    
    if(!allYBinsGood){
      std::cout << "MIXMACHINE::CheckMachinesMatchBins(): Y-Binning mismatch." << std::endl;
      std::cout << " Precision level: " << precision << std::endl;
      std::cout << " Main machine: " << m_mixMachineName << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << std::endl;

      for(unsigned int bIY = 0; bIY < checkBinsY.size(); ++bIY){
	std::cout << " " << bIY << ": " << m_binsYVect[bIY] << ", " << checkBinsY[bIY] << std::endl;
      }
      
      return false;
    }

  }
  else{
    Int_t nBinsXToCheck = machineToCheck->GetNBinsX();

    if(nBinsXToCheck != m_nBinsX){
      std::cout << "MIXMACHINE::CHECKMACHINEMATCHBINS() ERROR: nBinsX differs." << std::endl;
      std::cout << " Main Machine: " << m_mixMachineName << ", nBinsX: " << m_nBinsX << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << ", nBinsX: "  << nBinsXToCheck << std::endl;
      std::cout << "return false" << std::endl;
      return false;
    }

    std::vector<double> checkBinsX = machineToCheck->GetBinsX();
    bool allXBinsGood = true;
    for(unsigned int bIX = 0; bIX < checkBinsX.size(); ++bIX){
      if(TMath::Abs(checkBinsX[bIX] - m_binsXVect[bIX]) > precision){
        allXBinsGood = false;
        break;
      }
    }

    if(!allXBinsGood){
      std::cout << "MIXMACHINE::CheckMachinesMatchBins(): X-Binning mismatch." << std::endl;
      std::cout << " Precision level: " << precision << std::endl;
      std::cout << " Main machine: " << m_mixMachineName << std::endl;
      std::cout << " Machine-to-check: " << machineToCheckName << std::endl;

      for(unsigned int bIX = 0; bIX < checkBinsX.size(); ++bIX){
	std::cout << " " << bIX << ": " << m_binsXVect[bIX] << ", " << checkBinsX[bIX] << std::endl;
      }

      return false;
    }
  }

  return true;
}

TH1F* mixMachine::GetTH1FPtr(std::string histType)
{
  if(m_is2DUnfold){
    std::cout << "ERROR IN MIXMACHINE::GetTH1FPtr(): Unfold is 2-D; Calling GetTH1FPtr is invalid. return nullptr" << std::endl;
    return nullptr;
  }

  std::vector<std::string> mixNames = m_noneMixNames;
  if(m_mixMode == INCLUSIVE) mixNames = m_inclusiveMixNames;
  else if(m_mixMode == MULTI) mixNames = m_multiMixNames;

  int histPos = vectContainsStrPos(histType, &mixNames);
  if(histPos < 0){
    std::cout << "ERROR IN MIXMACHINE::GetTH1FPtr(): Requested hist \'" << histType << "\' is not found in MixMode \'" << m_mixMode << ". return nullptr" << std::endl;
    return nullptr;
  }
  
  return m_hists1D[histPos];
}

TH2F* mixMachine::GetTH2FPtr(std::string histType)
{
  if(!m_is2DUnfold){
    std::cout << "ERROR IN MIXMACHINE::GetTH2FPtr(): Unfold is 1-D; Calling GetTH2FPtr is invalid. return nullptr" << std::endl;
    return nullptr;
  }

  std::vector<std::string> mixNames = m_noneMixNames;
  if(m_mixMode == INCLUSIVE) mixNames = m_inclusiveMixNames;
  else if(m_mixMode == MULTI) mixNames = m_multiMixNames;

  int histPos = vectContainsStrPos(histType, &mixNames);

  if(histPos < 0){
    std::cout << "ERROR IN MIXMACHINE::GetTH2FPtr(): Requested hist \'" << histType << "\' is not found in MixMode \'" << m_mixMode << ". return nullptr" << std::endl;
    return nullptr;
  }
  
  return m_hists2D[histPos];
}

bool mixMachine::Add(mixMachine *machineToAdd, double precision)
{
  //Make sure the machines are compatible
  if(!(this->CheckMachinesMatchMode(machineToAdd))) return false;
  if(!(this->CheckMachinesMatchMC(machineToAdd))) return false;
  if(!(this->CheckMachinesMatch2D(machineToAdd))) return false;
  if(!(this->CheckMachinesMatchBins(machineToAdd, precision))) return false;

  if(m_is2DUnfold){
    std::vector<TH2F*> tempHists = machineToAdd->GetTH2F();
    for(unsigned int i = 0; i < tempHists.size(); ++i){
      m_hists2D[i]->Add(tempHists[i]);
    }
  }
  else{
    std::vector<TH1F*> tempHists = machineToAdd->GetTH1F();
    for(unsigned int i = 0; i < tempHists.size(); ++i){
      m_hists1D[i]->Add(tempHists[i]);
    }
  }

  return true;
}

bool mixMachine::Add(mixMachine* machineToAdd1, mixMachine* machineToAdd2, double precision)
{
  if(!this->Add(machineToAdd1, precision)) return false;
  if(!this->Add(machineToAdd2, precision)) return false;
  return true;
}

void mixMachine::ComputeSub()
{
  //  std::cout << "COMPUTING SUB FOR MIXMACHINE: \'" << m_mixMachineName << std::endl;

  if(m_mixMode == NONE){    
    //   std::cout << " MIXMODE IS NONE" << std::endl;
    int rawPos = vectContainsStrPos("RAW", &m_noneMixNames);
    int subPos = vectContainsStrPos("SUB", &m_noneMixNames);
    /*
    std::cout << " RAW POS: "<< rawPos << std::endl;
    std::cout << " SUB POS: "<< subPos << std::endl;
    */

    if(m_is2DUnfold) m_hists2D[subPos]->Add(m_hists2D[rawPos], 1.0);
    else m_hists1D[subPos]->Add(m_hists1D[rawPos], 1.0);
  }
  else if(m_mixMode == INCLUSIVE){
    //    std::cout << " MIXMODE IS INCLUSIVE" << std::endl;

    int rawPos = vectContainsStrPos("RAW", &m_inclusiveMixNames);
    int mixPos = vectContainsStrPos("MIX", &m_inclusiveMixNames);
    int subPos = vectContainsStrPos("SUB", &m_inclusiveMixNames);
    /*
    std::cout << " RAW POS: "<< rawPos << std::endl;
    std::cout << " MIX POS: "<< mixPos << std::endl;
    std::cout << " SUB POS: "<< subPos << std::endl;
    */
    if(m_is2DUnfold){
      m_hists2D[subPos]->Add(m_hists2D[rawPos], m_hists2D[mixPos], 1.0, -1.0);

      /*
      std::cout << "RAW HIST" << std::endl;
      for(Int_t bIY = m_hists2D[rawPos]->GetXaxis()->GetNbins()-1; bIY >= 0; --bIY){
	std::string outStr = "";
	
	for(Int_t bIX = 0; bIX < m_hists2D[rawPos]->GetXaxis()->GetNbins()+1; ++bIX){
	  std::string newStr = prettyString(m_hists2D[rawPos]->GetBinContent(bIX+1, bIY+1), 5, false);
	  while(newStr.size() > 7){newStr = newStr.replace(newStr.size()-1,1,"");}
	  outStr = outStr + newStr + " ";
	}
	std::cout << outStr << std::endl;
      }
      
      std::cout << "MIX HIST" << std::endl;
      for(Int_t bIY = m_hists2D[mixPos]->GetXaxis()->GetNbins()-1; bIY >= 0; --bIY){
	std::string outStr = "";
	
	for(Int_t bIX = 0; bIX < m_hists2D[mixPos]->GetXaxis()->GetNbins()+1; ++bIX){
	  std::string newStr = prettyString(m_hists2D[mixPos]->GetBinContent(bIX+1, bIY+1), 5, false);
	  while(newStr.size() > 7){newStr = newStr.replace(newStr.size()-1,1,"");}
	  outStr = outStr + newStr + " ";
	}
	std::cout << outStr << std::endl;
      }

      std::cout << "SUB HIST" << std::endl;
      for(Int_t bIY = m_hists2D[subPos]->GetXaxis()->GetNbins()-1; bIY >= 0; --bIY){
	std::string outStr = "";
	
	for(Int_t bIX = 0; bIX < m_hists2D[subPos]->GetXaxis()->GetNbins()+1; ++bIX){
	  std::string newStr = prettyString(m_hists2D[subPos]->GetBinContent(bIX+1, bIY+1), 5, false);
	  while(newStr.size() > 7){newStr = newStr.replace(newStr.size()-1,1,"");}
	  outStr = outStr + newStr + " ";
	}
	std::cout << outStr << std::endl;
      }

      std::cout << std::endl;
      */
    }
    else m_hists1D[subPos]->Add(m_hists1D[rawPos], m_hists1D[mixPos], 1.0, -1.0);
  }
  else if(m_mixMode == MULTI){
    std::cout << " MIXMODE IS MULTI" << std::endl;

    int rawPos = vectContainsStrPos("RAW", &m_multiMixNames);
    int mixPos = vectContainsStrPos("MIX", &m_multiMixNames);
    int mixCorrectionPos = vectContainsStrPos("MIXCORRECTION", &m_multiMixNames);
    int mixCorrectedPos = vectContainsStrPos("MIXCORRECTED", &m_multiMixNames);
    int subPos = vectContainsStrPos("SUB", &m_multiMixNames);
    std::cout << " RAW POS: "<< rawPos << std::endl;
    std::cout << " MIX POS: "<< mixPos << std::endl;
    std::cout << " MIXCORRECTION POS: "<< mixCorrectionPos << std::endl;
    std::cout << " MIXCORRECTED POS: "<< mixCorrectedPos << std::endl;
    std::cout << " SUB POS: "<< subPos << std::endl;

    if(m_is2DUnfold){
      m_hists2D[mixCorrectedPos]->Add(m_hists2D[mixPos], m_hists2D[mixCorrectionPos], 1.0, -1.0);
      m_hists2D[subPos]->Add(m_hists2D[rawPos], m_hists2D[mixCorrectedPos], 1.0, -1.0);
    }
    else{
      m_hists1D[mixCorrectedPos]->Add(m_hists1D[mixPos], m_hists1D[mixCorrectionPos], 1.0, -1.0);
      m_hists1D[subPos]->Add(m_hists1D[rawPos], m_hists1D[mixCorrectedPos], 1.0, -1.0);
    }
  }

  return;  
}


void mixMachine::WriteToFile(TFile* inFile_p)
{
  inFile_p->cd();
  
  if(m_is2DUnfold){
    for(unsigned int i = 0; i < m_hists2D.size(); ++i){
      m_hists2D[i]->Write("", TObject::kOverwrite);
    }    
  }
  else{
    for(unsigned int i = 0; i < m_hists1D.size(); ++i){
      m_hists1D[i]->Write("", TObject::kOverwrite);
    }    
  }

  return;
}

void mixMachine::WriteToDirectory(TDirectoryFile* inDir_p)
{
  inDir_p->cd();
  
  if(m_is2DUnfold){
    for(unsigned int i = 0; i < m_hists2D.size(); ++i){
      m_hists2D[i]->Write("", TObject::kOverwrite);
    }    
  }
  else{
    for(unsigned int i = 0; i < m_hists1D.size(); ++i){
      m_hists1D[i]->Write("", TObject::kOverwrite);
    }    
  }

  return;
}

void mixMachine::Print()
{
  std::cout << "mixMachine::Print():" << std::endl;
  std::cout << " mixMachineName: " << m_mixMachineName << std::endl;
  std::cout << " mixMode: " << m_mixMode << std::endl;
  std::cout << " is2DUnfold: " << m_is2DUnfold << std::endl;
  std::cout << " nBinsX: " << m_nBinsX << std::endl;
  std::cout << " binsX: ";
  for(Int_t bIX = 0; bIX < m_nBinsX; ++bIX){
    std::cout << m_binsX[bIX] << ", ";
  }
  std::cout << m_binsX[m_nBinsX] << "." << std::endl;
  std::cout << " titleX: " << m_titleX << std::endl;

  if(m_is2DUnfold){
    std::cout << " nBinsY: " << m_nBinsY << std::endl;
    std::cout << " binsY: ";
    for(Int_t bIY = 0; bIY < m_nBinsY; ++bIY){
      std::cout << m_binsY[bIY] << ", ";
    }
    std::cout << m_binsY[m_nBinsY] << "." << std::endl;
    std::cout << " titleY: " << m_titleY << std::endl;

    std::cout << std::endl;
    std::cout << " Histogram Names: " << std::endl;
    for(unsigned int i = 0; i < m_hists2D.size(); ++i){
      std::cout << "  " << i << "/" << m_hists2D.size() << ": " << m_hists2D[i]->GetName() << std::endl;
    }
  }
  else{
    std::cout << std::endl;
    std::cout << " Histogram Names: " << std::endl;
    for(unsigned int i = 0; i < m_hists1D.size(); ++i){
      std::cout << "  " << i << "/" << m_hists1D.size() << ": " << m_hists1D[i]->GetName() << std::endl;
    }
  }
  std::cout << "End mixMachine::Print()" << std::endl;

  return;
}

void mixMachine::Clean()
{
  m_isInit = false;
  m_mixMachineName = "";
  m_mixMode = NONE;

  m_env.Clear();

  m_nBinsX = 0;
  m_binsXVect.clear();
  m_nBinsY = 0;
  m_binsYVect.clear();
  m_titleX = "";
  m_titleY = "";

  for(unsigned int i = 0; i < m_hists1D.size(); ++i){
    delete m_hists1D[i];
  }  
  m_hists1D.clear();

  for(unsigned int i = 0; i < m_hists2D.size(); ++i){
    delete m_hists2D[i];
  }  
  m_hists2D.clear();

  return;
}
