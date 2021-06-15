//cpp dependencies
#include <iostream>

//Local dependencies
#include "include/mixMachine.h"

mixMachine::mixMachine(std::string inMixMachineName, mixMachine::mixMode inMixMode)
{
  if(!Init(inMixMachineName, inMixMode)) std::cout << "MIXMACHINE: Initialization failure. return" << std::endl;
  return;
}

bool mixMachine::Init(std::string inMixMachineName, mixMachine::mixMode inMixMode)
{
  m_mixMachineName = inMixMachineName;
  m_mixMode = inMixMode;

  m_isInit = true;
  return m_isInit;
}

void mixMachine::Clean()
{
  m_isInit = false;
  m_mixMachineName = "";
  m_mixMode = NONE;
}
