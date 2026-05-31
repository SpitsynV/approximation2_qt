#ifndef TASK31_H
#define TASK31_H

#include <vector>

// Вычисление всех производных в узлах сетки по методу "разделённые разности"
void computeAllDerivatives31(int nx, int ny,
                             const std::vector<double>& x,
                             const std::vector<double>& y,
                             const std::vector<double>& f_vals,
                             std::vector<double>& dx,
                             std::vector<double>& dy,
                             std::vector<double>& dxy);

// Вычисление значения интерполянта в точке
double GetValue31(double px, double py,
                  const std::vector<double>& x, int nx,
                  const std::vector<double>& y, int ny,
                  const std::vector<double>& f_vals,
                  const std::vector<double>& dx,
                  const std::vector<double>& dy,
                  const std::vector<double>& dxy);

#endif