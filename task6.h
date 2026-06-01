#pragma once

#include <vector>
#include <functional>

// Построение коэффициентов разложения по многочленам Чебышева (МНК).
// Параметры:
//   nx, ny  – количество узлов Чебышева по x и y
//   a, b    – границы по x
//   c, d    – границы по y
//   func    – функция f(x,y)
//   coeffs  – выходной массив коэффициентов (размер nx*ny)
void GetCoefficients6(int nx, int ny,
                      double a, double b, double c, double d,
                      std::function<double(double, double)> func,
                      std::vector<double>& coeffs);

// Вычисление значения разложения в произвольной точке (x,y).
double GetValue6(double x, double y,
                 double a, double b, double c, double d,
                 const std::vector<double>& coeffs,
                 int nx, int ny);
