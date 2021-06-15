#ifndef MIXMACHINE_H
#define MIXMACHINE_H

//cpp dependencies
#include <string>
#include <vector>

//ROOT dependencies
#include "TH1F.h"
#include "TH2F.h"

class mixMachine{
 public:
  enum mixMode{NONE=0,//Dummy mode - if on things should fail
	       INCLUSIVE= 1,
	       MULTI=2	  
  };

  mixMachine(){m_isInit = false;};
  ~mixMachine(){};

  mixMachine(std::string inMixMachineName, mixMachine::mixMode inMixMode);
  bool Init(std::string inMixMachineName, mixMachine::mixMode inMixMode);

  void Clean();

 private:
  bool m_isInit;
  std::string m_mixMachineName;
  mixMachine::mixMode m_mixMode;
  
};

#endif
