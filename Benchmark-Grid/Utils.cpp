#include "Utils.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <Grid/Grid.h>


double usecond_precise()
{
  using namespace std::chrono;
  auto nsecs = duration_cast<nanoseconds>(Grid::GridClock::now() - Grid::theProgramStart);
  return nsecs.count() * 1e-3;
}

void time_statistics::statistics(std::vector<double> v)
{
  double sum = std::accumulate(v.begin(), v.end(), 0.0);
  mean = sum / v.size();

  std::vector<double> diff(v.size());
  std::transform(v.begin(), v.end(), diff.begin(), [=](double x) { return x - mean; });
  double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  err = std::sqrt(sq_sum / (v.size() * (v.size() - 1)));

  auto result = std::minmax_element(v.begin(), v.end());
  min = *result.first;
  max = *result.second;
}
