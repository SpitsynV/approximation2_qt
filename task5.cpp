#include"task5.h"
// Реализация задачи : кусочная интерполяция кубическими многочленами Эрмита. Использует точные значения производных
#include <vector>
#include <cmath>
#include <algorithm>

// ----------------------------------------------------------------------
// Вычисление значения кубического интерполянта Эрмита на отрезке [x0,x1]
// в точке x. Известны: y0 = f(x0), y1 = f(x1), dy0 = f'(x0), dy1 = f'(x1).
// ----------------------------------------------------------------------
static double hermiteInterp1D(double x0, double x1, double y0, double y1, double dy0, double dy1, double x)
{
    double h = x1 - x0;
    double t = (x - x0) / h;          // локальная координата в [0,1]
    double t2 = t * t, t3 = t2 * t;

    // базисные функции Эрмита
    double H00 = 2*t3 - 3*t2 + 1;   
    double H10 = t3 - 2*t2 + t;     
    double H01 = -2*t3 + 3*t2;      
    double H11 = t3 - t2;           

    return y0 * H00 + dy0 * h * H10 + y1 * H01 + dy1 * h * H11;
}


// ----------------------------------------------------------------------
// Вычисление значения бикубического интерполянта в произвольной точке (px, py)
// с использованием последовательных одномерных интерполяций.
// Алгоритм:
//   Для каждого j (0..ny-1):
//     1. Находим отрезок по x, содержащий px.
//     2. Вычисляем val_j = интерполяция по x значений f и dx (в узлах этого отрезка).
//     3. Вычисляем der_y_j = интерполяция по x значений dy и dxy (в узлах этого отрезка).
//   Затем интерполируем по y: по точкам y_j, значениям val_j и производным der_y_j
//   получаем итоговое значение в py.
// Параметры:
//   px, py – координаты точки
//   x, y   – узлы сетки (размеры nx, ny)
//   f_vals – значения функции в узлах (плоский массив)
//   dx, dy, dxy – предвычисленные производные
// ----------------------------------------------------------------------
double GetValue5(double px, double py,
                  const std::vector<double>& x, int nx,
                  const std::vector<double>& y, int ny,
                  const std::vector<double>& f_vals,
                  const std::vector<double>& dx,
                  const std::vector<double>& dy,
                  const std::vector<double>& dxy)
{
    // 1. Находим отрезок по x, содержащий px (индекс левой границы i)
    int ix = std::lower_bound(x.begin(), x.end(), px) - x.begin();
    if (ix == 0) ix = 0;
    else if (ix == nx) ix = nx - 2;
    else ix = ix - 1;
    ix = std::max(0, std::min(nx-2, ix));
    double x0 = x[ix], x1 = x[ix+1];

    // 2. Для каждого y_j вычисляем интерполированные значения
    std::vector<double> F_j(ny), Fy_j(ny);
    for (int j = 0; j < ny; ++j) {
        // Извлекаем данные для отрезка [x0,x1] при данном y_j
        double y0 = f_vals[ix * ny + j];
        double y1 = f_vals[(ix+1) * ny + j];
        double dy0 = dx[ix * ny + j];     // ∂f/∂x в левом узле
        double dy1 = dx[(ix+1) * ny + j]; // ∂f/∂x в правом узле
        F_j[j] = hermiteInterp1D(x0, x1, y0, y1, dy0, dy1, px);

        // Аналогично интерполируем dy (производную по y) и dxy
        double z0 = dy[ix * ny + j];
        double z1 = dy[(ix+1) * ny + j];
        double dz0 = dxy[ix * ny + j];
        double dz1 = dxy[(ix+1) * ny + j];
        Fy_j[j] = hermiteInterp1D(x0, x1, z0, z1, dz0, dz1, px);
    }

    // 3. Интерполяция по y
    int iy = std::lower_bound(y.begin(), y.end(), py) - y.begin();
    if (iy == 0) iy = 0;
    else if (iy == ny) iy = ny - 2;
    else iy = iy - 1;
    iy = std::max(0, std::min(ny-2, iy));
    double y0 = y[iy], y1 = y[iy+1];

    double v0 = F_j[iy];
    double v1 = F_j[iy+1];
    double dv0 = Fy_j[iy];
    double dv1 = Fy_j[iy+1];

    return hermiteInterp1D(y0, y1, v0, v1, dv0, dv1, py);
}