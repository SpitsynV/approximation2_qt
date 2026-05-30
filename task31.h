#pragma once
#include <vector>

// Метод 31: кусочная интерполяция кубическими многочленами
//           с использованием разделённых разностей (§19.6 Богачев)
//           и естественными граничными условиями   (§19.8.3 Богачев)
//           для функций двух переменных            (§1.6.3 Глава II)
//
// Общая схема (§19.1): на каждом отрезке [x_i, x_{i+1}] (в 2D — ячейке)
//   P_i(x) = c_{1,i} + c_{2,i}*(x-x_i) + c_{3,i}*(x-x_i)^2 + c_{4,i}*(x-x_i)^3  (ф.2, стр.68)
//
// Параметры d_i (§19.6, ф.(9) стр.71):
//   d_i = sign(f(x_i;x_{i+1})) * min(|f(x_{i-1};x_i)|, |f(x_i;x_{i+1})|),
//         если sign(f(x_{i-1};x_i)) = sign(f(x_i;x_{i+1}))
//   d_i = 0, иначе  (i=2,...,n-1)
//
// Естественные ГУ (§19.8.3, стр.76): (Pf)''(x_1) = 0, (Pf)''(x_n) = 0:
//   2*d_1 + d_2   = 3*f(x_1;x_2)      =>  d_1 = (3*f(x_1;x_2) - d_2)   / 2
//   d_{n-1}+2*d_n = 3*f(x_{n-1};x_n)  =>  d_n = (3*f(x_{n-1};x_n) - d_{n-1}) / 2
//
// 2D-версия (§1.6.3): тензорное произведение — на ячейке [x_i,x_{i+1}]x[y_j,y_{j+1}]
//   P_{ij}(x,y) = sum_{k=0..3} sum_{l=0..3} gamma[k][l] * (x-x_i)^k * (y-y_j)^l
//   Gamma = A^{-1}(h_i) * F_{ij} * (A^{-1}(h_j))^T   (ф. стр.106)
//
// Хранение: coef[ (i*(ny-1)+j)*16 + k*4+l ] = gamma_{k+1,l+1}
//           i=0..nx-2, j=0..ny-2, k,l=0..3

// Построение коэффициентов кусочного многочлена
// f_vals[i*ny+j] = f(x[i], y[j])
void GetCoefficients31(int nx, int ny,
                       const std::vector<double> &x,
                       const std::vector<double> &y,
                       const std::vector<double> &f_vals,
                       std::vector<double> &coef);

// Вычисление значения в точке (px, py) — бинарный поиск + схема Горнера
double GetValue31(double px, double py,
                  const std::vector<double> &x, int nx,
                  const std::vector<double> &y, int ny,
                  const std::vector<double> &coef);