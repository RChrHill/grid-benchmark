#pragma once

#include <cxxabi.h>
#include <string>
#include <vector>


// NOTE: Grid::GridClock is just a typedef to
// `std::chrono::high_resolution_clock`, but `Grid::usecond` rounds to
// microseconds (no idea why, probably wasnt ever relevant before), so we need
// our own wrapper here.
double usecond_precise();

struct time_statistics
{
  double mean;
  double err;
  double min;
  double max;

  void statistics(std::vector<double> v);
};


template<typename T>
std::string getClassName()
{
  int status;
  char* name = abi::__cxa_demangle(typeid(T).name(),0,0,&status);
  std::string out{name};
  free(name);
  return out.substr(0, out.find("<"));
}
