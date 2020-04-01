//Author: Chris McGinn (2020.03.09)
//Contact at chmc7718@colorado.edu or cffionn on skype for bugs

//c+cpp
#include <fstream>
#include <iostream>
#include <map>
#include <string>

//Local
#include "include/checkMakeDir.h"

int grlToTex(std::string inGRLFileName)
{
  checkMakeDir check;
  if(!check.checkFileExt(inGRLFileName, ".xml")) return 1;

  std::ifstream inFile(inGRLFileName.c_str());
  std::string tempStr;

  std::map<std::string, std::string> runToLumiStr;

  std::cout << "\\begin{table}[h!]" << std::endl;
  std::cout << "\\fontsize{10}{10}\\selectfont" << std::endl;
  std::cout << "\\begin{center}" << std::endl;
  std::cout << "\\caption{Placeholder}" << std::endl;
  std::cout << "\\label{tab:placeholder}" << std::endl;
  std::cout << "\\hspace{-1.2cm}" << std::endl;
  std::cout << "\\begin{tabular}{ l l }" << std::endl;
  std::cout << "Run & Lumiblocks \\\\ \\hline" << std::endl;
  
  std::string currRunStr = "";
  std::string currLumiStr = "";
  std::string lineStr = "";
  while(std::getline(inFile, tempStr)){
    //    std::cout << tempStr << std::endl;

    if(tempStr.find("<Run") != std::string::npos){
      tempStr.replace(0, tempStr.find(">")+1, "");
      tempStr.replace(tempStr.rfind("<"), tempStr.size(), "");

      currRunStr = tempStr;

    }
    else if(currRunStr.size() != 0 && tempStr.find("<LB") != std::string::npos){
      tempStr.replace(0, tempStr.find("\"")+1, "");
      tempStr.replace(tempStr.rfind("\""), tempStr.size(), "");
      std::string firstNum = tempStr.substr(0, tempStr.find("\""));
      std::string secondNum = tempStr;
      while(secondNum.find("\"") != std::string::npos){
	secondNum.replace(0, secondNum.find("\"")+1, "");
      }
      while(firstNum.size() < 3){firstNum = "0" + firstNum;}
      while(secondNum.size() < 3){secondNum = "0" + secondNum;}
      currLumiStr = currLumiStr + firstNum + "-" + secondNum + ", ";
    }
    else if(tempStr.find("/LumiBlockCollection") != std::string::npos && currRunStr.size() != 0){
      currLumiStr.replace(currLumiStr.rfind(","), currLumiStr.size(), "");      
      if(lineStr.size() != 0) std::cout << lineStr << " \\\\" << std::endl;
      lineStr = currRunStr + " & " + currLumiStr;
      currRunStr = "";
      currLumiStr = "";
    }
  }  
  inFile.close();  

  std::cout << lineStr << std::endl;
  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\end{center}" << std::endl;
  std::cout << "\\end{table}" << std::endl;

  std::cout << "GRLTOTEX COMPLETE. return 0." << std::endl;
  return 0;
}

int main(int argc, char* argv[])
{
  if(argc != 2){
    std::cout << "Usage: ./bin/grlToTex.exe <inGRLFileName>" << std::endl;
    std::cout << "TO DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=1 #from command line" << std::endl;
    std::cout << "TO TURN OFF DEBUG:" << std::endl;
    std::cout << " export DOGLOBALDEBUGROOT=0 #from command line" << std::endl;
    std::cout << "return 1." << std::endl;
    return 1;
  }
 
  int retVal = 0;
  retVal += grlToTex(argv[1]);
  return retVal;
}
