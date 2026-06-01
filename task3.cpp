// task3.cpp --НАДО ДЛЯ БЕССЕЛЯ!
#include"task3.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ----------------------------------------------------------------------
// Одномерное вычисление производных d[i] для сетки x[0..n-1] и значений f[0..n-1]
// с использованием естественных граничных условий (P''=0 на концах).
// ----------------------------------------------------------------------
static void computeDerivatives1D(const double* x, const double* f, int n, double* d)
{
    if (n < 2) return;

    // разделённые разности первого порядка
    std::vector<double> dd(n - 1);
    for (int i = 0; i < n - 1; ++i)
        dd[i] = (f[i+1] - f[i]) / (x[i+1] - x[i]);

    if (n == 2) {
        d[0] = d[1] = dd[0];
        return;
    }

    // внутренние узлы:(метод Бесселя)
    for (int i = 1; i < n - 1; ++i)
        d[i] = ((x[i+1]-x[i])*d[i-1]+(x[i]-x[i-1])*d[i])/(x[i+1]-x[i-1]);

    // левая граница: 2*d[0] + d[1] = 3*dd[0]   → d[0] = (3*dd[0] - d[1])/2
    d[0] = (3.0 * dd[0] - d[1]) / 2.0;

    // правая граница: d[n-2] + 2*d[n-1] = 3*dd[n-2] → d[n-1] = (3*dd[n-2] - d[n-2])/2
    d[n-1] = (3.0 * dd[n-2] - d[n-2]) / 2.0;
}

// ----------------------------------------------------------------------
// Вычисление значения кубического интерполянта Эрмита на отрезке [x0,x1]
// в точке x. Известны: y0 = f(x0), y1 = f(x1), dy0 = f'(x0), dy1 = f'(x1).
// Используются формулы (3) на стр. 68.
// Возвращает интерполированное значение.
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
// Построение всех производных в узлах сетки (dx, dy, dxy).
// Вход: nx, ny, x, y, f_vals (плоский массив размера nx*ny).
// Выход: dx, dy, dxy – такие же плоские массивы.
// ----------------------------------------------------------------------
void computeAllDerivatives(int nx, int ny,
                                  const std::vector<double>& x,
                                  const std::vector<double>& y,
                                  const std::vector<double>& f_vals,
                                  std::vector<double>& dx,
                                  std::vector<double>& dy,
                                  std::vector<double>& dxy)
{
    /*надо ли? */
    dx.resize(nx * ny);
    dy.resize(nx * ny);
    dxy.resize(nx * ny);
    /*          */

    // 1. Производные по x (dx) для каждой строки j
    std::vector<double> f_col(nx);
    std::vector<double> d_col(nx);

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i)
            f_col[i] = f_vals[i * ny + j];
        computeDerivatives1D(x.data(), f_col.data(), nx, d_col.data());
        for (int i = 0; i < nx; ++i)
            dx[i * ny + j] = d_col[i];
    }

    // 2. Производные по y (dy) для каждого столбца i
    std::vector<double> f_row(ny);
    std::vector<double> d_row(ny);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j)
            f_row[j] = f_vals[i * ny + j];
        computeDerivatives1D(y.data(), f_row.data(), ny, d_row.data());
        for (int j = 0; j < ny; ++j)
            dy[i * ny + j] = d_row[j];
    }

    // 3. Смешанные производные dxy: применяем computeDerivatives1D к dy по x
    std::vector<double> dy_col(nx);
    std::vector<double> dxy_col(nx);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i)
            dy_col[i] = dy[i * ny + j];
        computeDerivatives1D(x.data(), dy_col.data(), nx, dxy_col.data());
        for (int i = 0; i < nx; ++i)
            dxy[i * ny + j] = dxy_col[i];
    }
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
double GetValue13(double px, double py,
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