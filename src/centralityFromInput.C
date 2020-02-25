//c+cpp
#include <fstream>

//Local
#include "include/centralityFromInput.h"
#include "include/plotUtilities.h"

centralityFromInput::centralityFromInput(std::string inTableFile){SetTable(inTableFile);}

void centralityFromInput::SetTable(std::string inTableFile)
{
  m_tableFileName = inTableFile;

  if(!m_check.checkFileExt(m_tableFileName, "txt")){
    std::cout << "CENTRALITYFROMINPUT: Given table \'" << m_tableFileName << "\' is invalid. return isInit=false" << std::endl;
    m_isInit = false;
    return;
  }

  std::ifstream inFile(m_tableFileName.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    while(tempStr.find(",") != std::string::npos){tempStr.replace(tempStr.find(","), 1, "");}
    if(tempStr.size() == 0) continue;

    m_centVals.push_back(std::stod(tempStr));

    if(m_centVals.size() > 1){
      if(m_centVals.size() == 2){
	if(m_centVals[m_centVals.size()-1] <= m_centVals[m_centVals.size()-2]) m_isDescending = true;
	else m_isDescending = false;
      }
      else{
	if(m_isDescending && m_centVals[m_centVals.size()-1] >= m_centVals[m_centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << m_tableFileName << "\' are not descending. return isInit=false" << std::endl;	
	  m_isInit = false;
	  return;
	}
	if(!m_isDescending && m_centVals[m_centVals.size()-1] <= m_centVals[m_centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << m_tableFileName << "\' are not ascending. return isInit=false" << std::endl;	
	  m_isInit = false;
	  return;
	}
      }
    }
  }  

  if(m_centVals.size() != 101){
    std::cout << "CENTRALITYFROMINPUT: Values in table \'" << m_tableFileName << "\' are not N=101. return isInit=false" << std::endl;	
    m_isInit = false;
    return;
  }
  
  m_isInit = true;
  
  return;
}

double centralityFromInput::GetCent(double inVal)
{
  double outVal = -1;

  if(!m_isInit) std::cout << "CENTRALITYFROMINPUT: Initialization failed. GetCent call will return -1" << std::endl;
  else{
    for(unsigned int cI = 0; cI < m_centVals.size()-1; ++cI){
      if(m_isDescending){
	if(inVal < m_centVals[cI] && inVal >= m_centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
      else{
	if(inVal > m_centVals[cI] && inVal <= m_centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
    }
  }

  return outVal;
}

void centralityFromInput::PrintTableTex()
{
  if(!m_isInit){
    std::cout << "CENTRALITYFROMINPUT ERROR - PrintTableTex(): Table is not initialized. return" << std::endl;
    return;
  }

  std::string tempFileName = m_tableFileName;
  unsigned int pos = 0;
  while(pos < tempFileName.size()){
    bool sizeInc = false;

    if(tempFileName.substr(pos,1).find("_") != std::string::npos){
      if(pos != 0){
	if(tempFileName.substr(pos-1, 1).find("\\") == std::string::npos){
	  tempFileName = tempFileName.substr(0, pos) + "\\" + tempFileName.substr(pos, tempFileName.size());
	  sizeInc = true;
	}	
      }
      else{
	tempFileName = "\\" + tempFileName;
	sizeInc = true;	     
      }
    }

    if(!sizeInc) ++pos;
  }				 
  std::cout << "TEMPFILENAME: " << tempFileName << std::endl;

    
  std::cout << "\\begin{table}[h!]" << std::endl;
  std::cout << "\\fontsize{8}{8}\\selectfont" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\caption{FCal $\\Sigma$E$_{\\mathrm{T}}$ to centrality as defined in file '" << tempFileName << "'}" << std::endl;

  std::cout << "\\label{tab:fcalToCent}" << std::endl;
  std::cout << "\\hspace{-1.2cm}" << std::endl;
  std::cout << "\\begin{tabular}{ l l | l l | l l }" << std::endl;
  std::cout << "Centrality (\\%) & FCal $\\Sigma$E$_{\\mathrm{T}}$ [GeV] & Centrality (\\%) & FCal $\\Sigma$E$_{\\mathrm{T}}$ [GeV] & Centrality (\\%) & FCal $\\Sigma$E$_{\\mathrm{T}}$ [GeV] \\\\ \\hline" << std::endl;

  int offsetVal = 99;
  if(m_isDescending) offsetVal = 0;

  const int nCent = 33;
  for(int cI = 0; cI < nCent; ++cI){
    double firstVal = m_centVals[cI];
    if(cI == 0) firstVal = 0;
    
    std::string str1 = std::to_string(offsetVal-cI) + "-" + std::to_string(offsetVal-cI+1) + "\\% & " + prettyString(firstVal, 3, false) + "-" + prettyString(m_centVals[cI+1], 3, false);
    std::string str2 = std::to_string(offsetVal - (nCent+cI)) + "-" + std::to_string(offsetVal - (nCent+cI)+1) + "\\% & " + prettyString(m_centVals[nCent+cI], 3, false) + "-" + prettyString(m_centVals[nCent+cI+1], 3, false);
    std::string str3 = std::to_string(offsetVal - (nCent*2+cI)) + "-" + std::to_string(offsetVal - (nCent*2+cI)+1) + "\\% & " + prettyString(m_centVals[nCent*2+cI], 3, false) + "-" + prettyString(m_centVals[nCent*2+cI+1], 3, false);

    std::cout << str1 << " & " << str2 << " & " << str3 << " \\\\" << std::endl;
  }
  std::cout << " & & & & 0-1\\%" << " & " << prettyString(m_centVals[m_centVals.size()-2], 3, false) + "-" + prettyString(m_centVals[m_centVals.size()-1], 3, false) << std::endl;

  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\end{table}" << std::endl;
  
  return;
}
