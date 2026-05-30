#pragma once
#include <vector>
double GetValue13(double px, double py,
                  const std::vector<double>& x, int nx,
                  const std::vector<double>& y, int ny,
                  const std::vector<double>& f_vals,
                  const std::vector<double>& dx,
                  const std::vector<double>& dy,
                  const std::vector<double>& dxy);

void computeAllDerivatives(int nx, int ny,
                                  const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& f_vals,
                                  std::vector<double>& dx,
                                  std::vector<double>& dy,
                                  std::vector<double>& dxy);
