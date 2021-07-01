//cpp dependencies
#include <iostream>

//ROOT dependencies
#include "THashList.h"

//Local dependencies
#include "include/envUtil.h"
#include "include/mixMachine.h"

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
  std::vector<Float_t> tempBins = strToVectF(m_env.GetValue("BINSX", ""));
  for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
    m_binsX[tI] = tempBins[tI];
  }
  m_titleX = m_env.GetValue("TITLEX", "");  

  if(m_is2DUnfold){
    if(!checkEnvForParams(&m_env, reqParams2DUnfold)){
      Clean();
      return false;
    }

    m_nBinsY = m_env.GetValue("NBINSY", -1);
    std::vector<Float_t> tempBins = strToVectF(m_env.GetValue("BINSY", ""));
    for(unsigned int tI = 0; tI < tempBins.size(); ++tI){
      m_binsY[tI] = tempBins[tI];
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

  if(m_is2DUnfold) m_hists2D[histPos]->Fill(fillX, fillY, fillWeight);
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

void mixMachine::ComputeSub()
{
  if(m_mixMode == NONE){
    int rawPos = vectContainsStrPos("RAW", &m_noneMixNames);
    int subPos = vectContainsStrPos("SUB", &m_noneMixNames);

    if(m_is2DUnfold) m_hists2D[subPos]->Add(m_hists2D[rawPos], 1.0);
    else m_hists1D[subPos]->Add(m_hists1D[rawPos], 1.0);
  }
  else if(m_mixMode == INCLUSIVE){
    int rawPos = vectContainsStrPos("RAW", &m_inclusiveMixNames);
    int mixPos = vectContainsStrPos("MIX", &m_inclusiveMixNames);
    int subPos = vectContainsStrPos("SUB", &m_inclusiveMixNames);

    if(m_is2DUnfold) m_hists2D[subPos]->Add(m_hists2D[rawPos], m_hists2D[mixPos], -1.0);
    else m_hists1D[subPos]->Add(m_hists1D[rawPos], m_hists1D[mixPos], -1.0);
  }
  else if(m_mixMode == MULTI){
    int rawPos = vectContainsStrPos("RAW", &m_multiMixNames);
    int mixPos = vectContainsStrPos("MIX", &m_multiMixNames);
    int mixCorrectionPos = vectContainsStrPos("MIXCORRECTION", &m_multiMixNames);
    int mixCorrectedPos = vectContainsStrPos("MIXCORRECTED", &m_multiMixNames);
    int subPos = vectContainsStrPos("SUB", &m_multiMixNames);

    if(m_is2DUnfold){
      m_hists2D[mixCorrectedPos]->Add(m_hists2D[mixPos], m_hists2D[mixCorrectionPos], -1.0);
      m_hists2D[subPos]->Add(m_hists2D[rawPos], m_hists2D[mixCorrectedPos], -1.0);
    }
    else{
      m_hists1D[mixCorrectedPos]->Add(m_hists1D[mixPos], m_hists1D[mixCorrectionPos], -1.0);
      m_hists1D[subPos]->Add(m_hists1D[rawPos], m_hists1D[mixCorrectedPos], -1.0);
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
  m_nBinsY = 0;
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
