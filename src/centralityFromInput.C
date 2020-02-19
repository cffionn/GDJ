//Local
#include "include/centralityFromInput.h"

centralityFromInput::centralityFromInput(std::string inTableFile){SetTable(inTableFile);}

void centralityFromInput::SetTable(std::string inTableFile)
{
  if(!check.checkFileExt(inTableFile, "txt")){
    std::cout << "CENTRALITYFROMINPUT: Given table \'" << inTableFile << "\' is invalid. return isInit=false" << std::endl;
    isInit = false;
    return;
  }

  std::ifstream inFile(inTableFile.c_str());
  std::string tempStr;
  while(std::getline(inFile, tempStr)){
    while(tempStr.find(",") != std::string::npos){tempStr.replace(tempStr.find(","), 1, "");}
    if(tempStr.size() == 0) continue;

    centVals.push_back(std::stod(tempStr));

    if(centVals.size() > 1){
      if(centVals.size() == 2){
	if(centVals[centVals.size()-1] <= centVals[centVals.size()-2]) isDescending = true;
	else isDescending = false;
      }
      else{
	if(isDescending && centVals[centVals.size()-1] >= centVals[centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not descending. return isInit=false" << std::endl;	
	  isInit = false;
	  return;
	}
	if(!isDescending && centVals[centVals.size()-1] <= centVals[centVals.size()-2]){
	  std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not ascending. return isInit=false" << std::endl;	
	  isInit = false;
	  return;
	}
      }
    }
  }  

  if(centVals.size() != 101){
    std::cout << "CENTRALITYFROMINPUT: Values in table \'" << inTableFile << "\' are not N=101. return isInit=false" << std::endl;	
    isInit = false;
    return;
  }
  
  isInit = true;
  
  return;
}

double centralityFromInput::GetCent(double inVal)
{
  double outVal = -1;

  if(!isInit) std::cout << "CENTRALITYFROMINPUT: Initialization failed. GetCent call will return -1" << std::endl;
  else{
    for(unsigned int cI = 0; cI < centVals.size()-1; ++cI){
      if(isDescending){
	if(inVal < centVals[cI] && inVal >= centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
      else{
	if(inVal > centVals[cI] && inVal <= centVals[cI+1]){
	  outVal = 99-cI;
	  break;
	}
      }
    }
  }

  return outVal;
}
