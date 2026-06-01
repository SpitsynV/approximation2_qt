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
// Вычисление значения в произвольной точке (px, py)
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
    // 1. Поиск отрезка по x
    int ix = std::lower_bound(x.begin(), x.end(), px) - x.begin();
    if (ix == 0) ix = 0;
    else if (ix == nx) ix = nx - 2;
    else ix = ix - 1;
    ix = std::max(0, std::min(nx - 2, ix));
    double x0 = x[ix], x1 = x[ix + 1];

    // 2. Поиск отрезка по y 
    int iy = std::lower_bound(y.begin(), y.end(), py) - y.begin();
    if (iy == 0) iy = 0;
    else if (iy == ny) iy = ny - 2;
    else iy = iy - 1;
    iy = std::max(0, std::min(ny - 2, iy));
    double y0 = y[iy], y1 = y[iy + 1];

    // 3. Интерполяция по x только для нужных j
    int idx00 = ix * ny + iy;         
    int idx10 = (ix + 1) * ny + iy;   

    double v0 = hermiteInterp1D(x0, x1,
        f_vals[idx00], f_vals[idx10],
        dx[idx00], dx[idx10], px);

    double dv0 = hermiteInterp1D(x0, x1,
        dy[idx00], dy[idx10],
        dxy[idx00], dxy[idx10], px);

    int idx01 = ix * ny + iy + 1;
    int idx11 = (ix + 1) * ny + iy + 1;

    double v1 = hermiteInterp1D(x0, x1,
        f_vals[idx01], f_vals[idx11],
        dx[idx01], dx[idx11], px);

    double dv1 = hermiteInterp1D(x0, x1,
        dy[idx01], dy[idx11],
        dxy[idx01], dxy[idx11], px);

    // 4. Финальная интерполяция по y
    return hermiteInterp1D(y0, y1, v0, v1, dv0, dv1, py);
}