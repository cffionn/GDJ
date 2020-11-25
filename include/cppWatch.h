#ifndef CPPWATCH_H
#define CPPWATCH_H

#include <ctime>

class cppWatch{
 public:
  cppWatch();
  ~cppWatch(){}

  double totalIntCPU;
  double currentIntCPU;
  double totalIntWall;
  double currentIntWall;
  std::clock_t c_start;
  std::time_t t_start;

  void start(){c_start = std::clock(); t_start = std::time(NULL); return;};
  void stop();
  double totalCPU(){return totalIntCPU;}
  double totalWall(){return totalIntWall;}
  double currentCPU(){return currentIntCPU;}
  double currentWall(){return currentIntWall;}
  void clear(){totalIntCPU = 0; currentIntCPU = 0; totalIntWall = 0; currentIntWall = 0; return;}
};

cppWatch::cppWatch(){totalIntCPU = 0; currentIntCPU = 0; totalIntWall = 0; currentIntWall = 0; return;}

void cppWatch::stop()
{
  std::clock_t c_end = std::clock();
  std::time_t t_end = std::time(NULL);

  totalIntWall += t_end  - t_start;
  currentIntWall = t_end  - t_start;
  totalIntCPU += c_end - c_start; 
  currentIntCPU = c_end - c_start;
  return;
}

#endif
