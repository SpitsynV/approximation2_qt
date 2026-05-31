#pragma once
#include <vector>
//ЭРМИТ 2D
//получаем ТОЧНЫЕ значения всех произодных
double GetValue5(double px, double py,
                  const std::vector<double>& x, int nx,
                  const std::vector<double>& y, int ny,
                  const std::vector<double>& f_vals,
                  const std::vector<double>& dx,
                  const std::vector<double>& dy,
                  const std::vector<double>& dxy);
