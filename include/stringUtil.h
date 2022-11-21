#ifndef STRINGUTIL_H
#define STRINGUTIL_H

//c+cpp
#include <ctime>
#include <iostream>
#include <map>
#include <string>
#include <vector>

//ROOT
#include "TEnv.h"
#include "TH1F.h"

inline bool isStrSame(std::string inStr1, std::string inStr2)
{
  if(inStr1.size() != inStr2.size()) return false;
  if(inStr1.find(inStr2) == std::string::npos) return false;
  return true;
}

inline std::string removeAllWhiteSpace(std::string inStr)
{
  while(inStr.find(" ") != std::string::npos){
    inStr.replace(inStr.find(" "), 1, "");
  }

  return inStr;
}

inline std::string returnAllCapsString(std::string inStr)
{
  const std::string lowStr = "abcdefghijklmnopqrstuvwxyz";
  const std::string hiStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  for(unsigned int lowIter = 0; lowIter < lowStr.size(); ++lowIter){
    while(inStr.find(lowStr.substr(lowIter, 1)) != std::string::npos){
      inStr.replace(inStr.find(lowStr.substr(lowIter, 1)), 1, hiStr.substr(lowIter, 1));
    }
  }

  return inStr;
}

inline std::string returnAllLowercaseString(std::string inStr)
{
  const std::string lowStr = "abcdefghijklmnopqrstuvwxyz";
  const std::string hiStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

  for(unsigned int hiIter = 0; hiIter < hiStr.size(); ++hiIter){
    while(inStr.find(hiStr.substr(hiIter, 1)) != std::string::npos){
      inStr.replace(inStr.find(hiStr.substr(hiIter, 1)), 1, lowStr.substr(hiIter, 1));
    }
  }

  return inStr;
}

inline bool isStrFromCharSet(const std::string inStr, const std::string charSet)
{
  for(unsigned int iter = 0; iter < inStr.size(); ++iter){
    if(charSet.find(inStr.substr(iter, 1)) == std::string::npos){
      return false;
    }
  }

  return true;
}

inline bool isStrAllAlpha(std::string inStr){return isStrFromCharSet(returnAllCapsString(inStr), "ABCDEFGHIJKLMNOPQRSTUVWXYZ");}
inline bool isStrInt(std::string inStr){return isStrFromCharSet(inStr, "-0123456789");}
inline bool isStrFloatOrDouble(std::string inStr){return isStrFromCharSet(inStr, ".-0123456789");}

inline bool isStrTrueOrFalse(std::string inStr)
{
  inStr = returnAllCapsString(inStr);
  if(!isStrAllAlpha(inStr)) return false;

  if(inStr.size() == 4 && inStr.find("TRUE") != std::string::npos) return true;
  if(inStr.size() == 5 && inStr.find("FALSE") != std::string::npos) return true;

  return false;
}

inline bool strToTrueOrFalse(std::string inStr)
{
  inStr = returnAllCapsString(inStr);

  if(inStr.size() == 4 && inStr.find("TRUE") != std::string::npos) return true;
  else if(inStr.size() == 5 && inStr.find("FALSE") != std::string::npos) return false;

  std::cout << "Call to strToTrueOrFalse is invalid; \'" << inStr << "\' neither true or false str. return false but you really ought to fix this" << std::endl;
  return false;
}

inline int getRVal(const std::string inStr)
{
  int rVal = -1;
  if(inStr.find("ak1PF") != std::string::npos) rVal = 1;
  else if(inStr.find("ak2PF") != std::string::npos) rVal = 2;
  else if(inStr.find("ak3PF") != std::string::npos) rVal = 3;
  else if(inStr.find("ak4PF") != std::string::npos) rVal = 4;
  else if(inStr.find("ak5PF") != std::string::npos) rVal = 5;
  else if(inStr.find("ak6PF") != std::string::npos) rVal = 6;
  else if(inStr.find("ak8PF") != std::string::npos) rVal = 8;
  else if(inStr.find("ak10PF") != std::string::npos) rVal = 10;
  else if(inStr.find("akCs1P") != std::string::npos) rVal = 1;
  else if(inStr.find("akCs2P") != std::string::npos) rVal = 2;
  else if(inStr.find("akCs3P") != std::string::npos) rVal = 3;
  else if(inStr.find("akCs4P") != std::string::npos) rVal = 4;
  else if(inStr.find("akCs5P") != std::string::npos) rVal = 5;
  else if(inStr.find("akCs6P") != std::string::npos) rVal = 6;
  else if(inStr.find("akCs8P") != std::string::npos) rVal = 8;
  else if(inStr.find("akCs10P") != std::string::npos) rVal = 10;
  else if(inStr.find("akPu3P") != std::string::npos) rVal = 3;
  else if(inStr.find("akPu4P") != std::string::npos) rVal = 4;

  return rVal;
}

inline std::string getRValStr(const std::string inStr)
{
  std::string rVal = "999";
  if(inStr.find("ak1PF") != std::string::npos) rVal = "0.1";
  else if(inStr.find("ak2PF") != std::string::npos) rVal = "0.2";
  else if(inStr.find("ak3PF") != std::string::npos) rVal = "0.3";
  else if(inStr.find("ak4PF") != std::string::npos) rVal = "0.4";
  else if(inStr.find("ak5PF") != std::string::npos) rVal = "0.5";
  else if(inStr.find("ak6PF") != std::string::npos) rVal = "0.6";
  else if(inStr.find("ak8PF") != std::string::npos) rVal = "0.8";
  else if(inStr.find("ak10PF") != std::string::npos) rVal = "1.0";
  else if(inStr.find("akCs1P") != std::string::npos) rVal = "0.1";
  else if(inStr.find("akCs2P") != std::string::npos) rVal = "0.2";
  else if(inStr.find("akCs3P") != std::string::npos) rVal = "0.3";
  else if(inStr.find("akCs4P") != std::string::npos) rVal = "0.4";
  else if(inStr.find("akCs5P") != std::string::npos) rVal = "0.5";
  else if(inStr.find("akCs6P") != std::string::npos) rVal = "0.6";
  else if(inStr.find("akCs8P") != std::string::npos) rVal = "0.8";
  else if(inStr.find("akCs10P") != std::string::npos) rVal = "1.0";
  else if(inStr.find("akPu3P") != std::string::npos) rVal = "0.3";
  else if(inStr.find("akPu4P") != std::string::npos) rVal = "0.4";

  return rVal;
}

inline std::vector<std::string> commaSepStringToVect(std::string inStr)
{
  std::vector<std::string> retVect;

  inStr = inStr + ",";
  while(inStr.find(",,") != std::string::npos){inStr.replace(inStr.find(",,"), 2, ",");}
  if(inStr.size() != 0){
    while(inStr.substr(0,1).find(",") != std::string::npos){
      inStr.replace(0, 1, "");
      if(inStr.size() == 0) break;
    }
  }

  while(inStr.find(",") != std::string::npos){
    retVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  
  return retVect;
}

inline bool vectContainsStr(std::string inStr, std::vector<std::string>* inVect)
{
  bool isInVect = false;
  if(inStr.size() != 0){
    for(unsigned int vI = 0; vI < inVect->size(); ++vI){
      if(isStrSame(inStr, (*inVect)[vI])){
	isInVect = true;
	break;
      }
    }
  }
  return isInVect;
}

inline int vectContainsStrPos(std::string inStr, std::vector<std::string>* inVect)
{
  int vectPos = -1;
  if(inStr.size() != 0){
    for(unsigned int vI = 0; vI < inVect->size(); ++vI){
      if(isStrSame(inStr, (*inVect)[vI])){
	vectPos = vI;
	break;
      }
    }
  }
  return vectPos;
}

inline bool vectContainsInt(int inInt, std::vector<int>* inVect)
{
  bool isInVect = false;
  for(unsigned int vI = 0; vI < inVect->size(); ++vI){
    if(inInt == ((*inVect)[vI])){
      isInVect = true;
      break;
    }
  }
  return isInVect;
}

inline bool vectContainsULL(unsigned long long inInt, std::vector<unsigned long long>* inVect)
{
  bool isInVect = false;
  for(unsigned int vI = 0; vI < inVect->size(); ++vI){
    if(inInt == ((*inVect)[vI])){
      isInVect = true;
      break;
    }
  }
  return isInVect;
}

inline std::vector<std::string> strToVect(std::string inStr)
{
  std::vector<std::string> retVect;
  while(inStr.find(",") != std::string::npos){
    retVect.push_back(inStr.substr(0, inStr.find(",")));
    inStr.replace(0, inStr.find(",")+1, "");
  }
  if(inStr.size() != 0) retVect.push_back(inStr);

  return retVect;
}

inline std::vector<float> strToVectF(std::string inStr)
{
  std::vector<float> retVect;
  std::vector<std::string> strVect = strToVect(inStr);
  for(auto const & iter : strVect){
    retVect.push_back(std::stod(iter));
  }
  return retVect;
}

inline std::vector<int> strToVectI(std::string inStr)
{
  std::vector<int> retVect;
  std::vector<std::string> strVect = strToVect(inStr);
  for(auto const & iter : strVect){
    retVect.push_back(std::stoi(iter));
  }
  return retVect;
}

inline std::vector<double> strToVectD(std::string inStr)
{
  std::vector<double> retVect;
  std::vector<std::string> strVect = strToVect(inStr);
  for(auto const & iter : strVect){
    retVect.push_back(std::stod(iter));
  }
  return retVect;
}

inline std::string vectToStr(std::vector<std::string> inVect, std::string sepChar)
{
  std::string retStr = "";
  for(unsigned int vI = 0; vI < inVect.size(); ++vI){
    retStr = retStr + inVect[vI] + sepChar;
  }
  return retStr;
}

inline std::string vectToStrComma(std::vector<std::string> inVect)
{
  return vectToStr(inVect, ",");
}

inline std::string getDateStr()
{
  std::time_t now = time(0);
  std::tm *ltm = std::localtime(&now);
  std::string dateStr = std::to_string(1+ltm->tm_mon);
  if(dateStr.size() == 1) dateStr = "0" + dateStr;
  std::string dayStr = std::to_string(ltm->tm_mday);
  if(dayStr.size() == 1) dayStr = "0" + dayStr;
  dateStr = std::to_string(1900 + ltm->tm_year) + dateStr + dayStr;

  return dateStr;
}

inline std::string rootFileNameProc(std::string inFileName, std::vector<std::string> modStrs)
{
  const std::string ext = ".root";
  const std::string defaultName = "stringutil_rootfilenameproc_default";
  std::string outFileName = inFileName;

  if(outFileName.size() == 0){
    std::cout << "STRINGUTIL ROOTFILENAMEPROC: Given inFileName is size 0, call defult name \'" << defaultName << "\'." << std::endl;
    outFileName = defaultName;
  }

  if(outFileName.size() >= ext.size()){
    if(isStrSame(outFileName, ext)){
      std::cout << "STRINGUTIL ROOTFILENAMEPROC: Given inFileName is size 0, call defult name \'" << defaultName << "\'." << std::endl;
      outFileName = defaultName;
    }
    else if(isStrSame(outFileName.substr(outFileName.size()-ext.size(), ext.size()), ext)) outFileName.replace(outFileName.size()-ext.size(), ext.size(), "");
  }

  for(unsigned int sI = 0; sI < modStrs.size(); ++sI){
    if(outFileName.size() < modStrs[sI].size() + 1) outFileName = outFileName + "_" + modStrs[sI];
    else if(!isStrSame(modStrs[sI], outFileName.substr(outFileName.size() - modStrs[sI].size(), modStrs[sI].size()))) outFileName = outFileName + "_" + modStrs[sI];
  }
  outFileName = outFileName + ".root";
  return outFileName;
}

inline std::string strLowerToUpper(std::string inStr)
{
  const std::string upper = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  const std::string lower = "abcdefghijklmnopqrstuvwxyz";

  for(unsigned int i = 0; i < inStr.size(); ++i){
    if(lower.find(inStr.substr(i, 1)) != std::string::npos){
      
      inStr.replace(i,1,upper.substr(lower.find(inStr.substr(i, 1)),1));
    }
  }
  return inStr;
}


template <typename T>
inline std::vector<std::string> getLabels(TEnv* plotConfig_p, T* histForLabels_p, std::map<std::string, std::string>* labelMap, std::vector<std::string>* labelsForSaveStr = nullptr)
{
  std::vector<std::string> labelVect;
  Int_t nGlobalLabels = 0;
  for(Int_t gI = 0; gI < 10; ++gI){
    std::string tempStr = plotConfig_p->GetValue(("MIXEDEVTPLOT.GLOBALLABEL." + std::to_string(gI)).c_str(), "");
    if(tempStr.size() == 0) break;
    while(tempStr.find("\"") != std::string::npos){tempStr.replace(tempStr.find("\""), 1, "");}

    labelVect.push_back(tempStr);
    ++nGlobalLabels;
  }

  int jetR = plotConfig_p->GetValue("JETR", 100);
  std::string jetRStr = "anti-k_{t} #it{R}=";
  if(jetR < 10) jetRStr = jetRStr + "0." + std::to_string(jetR) + " jets";
  //  else{
  //    jetRStr = jetRStr + prettyString(((double)jetR)/10., 1, false) + " jets";
  //  }
  labelVect.push_back(jetRStr);

  
  std::string histName = histForLabels_p->GetName();
  std::string labelStr = histName;
  if(labelStr.find("_") != std::string::npos) labelStr.replace(0, labelStr.find("_")+1, "");
  if(labelStr.rfind("_h") != std::string::npos) labelStr.replace(labelStr.rfind("_h"), 2, "");

  while(labelStr.find("_") != std::string::npos){
    std::string centStr = labelStr.substr(0, labelStr.find("_"));
    labelStr.replace(0, labelStr.find("_")+1, "");

    if(histName.find("DPhiVCent") != std::string::npos || histName.find("DPhiJJGVCent") != std::string::npos){
      if(labelStr.find("DPhi") == std::string::npos) continue;
    }

    if(centStr.find("Alt") == std::string::npos) labelVect.push_back(centStr);
    if(isStrSame("RAW_h", labelStr)) break;
  }
  if(labelStr.size() != 0) labelVect.push_back(labelStr);


  if(labelsForSaveStr != nullptr){
    for(unsigned int pI = nGlobalLabels; pI < labelVect.size(); ++pI){
      if(labelVect[pI].find("anti-k") != std::string::npos) continue;
      if(labelVect[pI].find("ATLAS") != std::string::npos) continue;

      labelsForSaveStr->push_back(labelVect[pI]);
    }
  }

  for(unsigned int pI = 0; pI < labelVect.size(); ++pI){
    if(labelMap->count(labelVect[pI]) != 0){
      labelVect[pI] = (*labelMap)[labelVect[pI]];
    }
  }


  std::string jtPtLowStr = plotConfig_p->GetValue("JTPTBINSLOW", "");
  std::string jtPtHighStr = plotConfig_p->GetValue("JTPTBINSHIGH", "");
  std::string jtPtLowRecoStr = plotConfig_p->GetValue("JTPTBINSLOWRECO", "");
  std::string jtPtHighRecoStr = plotConfig_p->GetValue("JTPTBINSHIGHRECO", "");
  std::string jtEtaLowStr = plotConfig_p->GetValue("JTETABINSLOW", "");
  std::string jtEtaHighStr = plotConfig_p->GetValue("JTETABINSHIGH", "");
  std::string mixJJDR = plotConfig_p->GetValue("MIXJETEXCLUSIONDR", "");
  //Add the three relevant jet cuts
  labelVect.push_back(jtPtLowStr + "<p_{T,Jet}^{Truth}<" + jtPtHighStr);
  labelVect.push_back(jtPtLowRecoStr + "<p_{T,Jet}^{Reco.}<" + jtPtHighRecoStr);
  labelVect.push_back(jtEtaLowStr + "<#eta_{Jet}<" + jtEtaHighStr);
  labelVect.push_back("#DeltaR_{JJ}>" + mixJJDR);

    for(unsigned int pI = 0; pI < labelVect.size(); ++pI){
    if(labelVect[pI].find("Barrel") != std::string::npos || labelVect[pI].find("EC") != std::string::npos) labelVect[pI] = labelVect[pI] + " #gamma";
  }

  if(histName.find("XJJVCent") != std::string::npos || histName.find("AJJVCent") != std::string::npos){
    std::string multiDPhiJJStr = "|#Delta#phi_{#vec{JJ},#gamma}| > " + std::string(plotConfig_p->GetValue("GAMMAMULTIJTDPHI", ""));
    //    std::cout << "MULTIJTDPHI 1: " << multiDPhiJJStr << std::endl;                                
    if(multiDPhiJJStr.find("pi") != std::string::npos) multiDPhiJJStr.replace(multiDPhiJJStr.find("pi"), 2, "#pi");
    //    std::cout << "MULTIJTDPHI 2: " << multiDPhiJJStr << std::endl;                                
    labelVect.push_back(multiDPhiJJStr);
  }

  /*
  const bool isMC = plotConfig_p->GetValue("ISMC", 0);
  if(isMC){
    std::string assocGenMinPtStr = plotConfig_p->GetValue("ASSOCGENMINPT", "");
    if(assocGenMinPtStr.size() != 0) labelVect.push_back("p_{T,Jet}^{Truth Match} > " + assocGenMinPtStr);
  }
  */

  return labelVect;

}

inline bool strReplace(std::string* inStr, std::string strToRep, std::string repStr)
{
  bool goodReplace = true;
  if(inStr->find(strToRep) != std::string::npos) inStr->replace(inStr->find(strToRep), strToRep.size(), repStr);
  else{
    std::cout << "STRREPLACE ERROR: String '" << *inStr << "' does not contain '" << strToRep << "' requested for replacement. return false" << std::endl;
    goodReplace = false;
  }

  return goodReplace;
}

#endif
