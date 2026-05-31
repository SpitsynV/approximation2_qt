#include "task31.h" //по методу раздел разностей -ОК
#include <vector>
#include <cmath>
#include <algorithm>

// ----------------------------------------------------------------------
// Одномерное вычисление производных по правилу "разделённые разности"
// с естественными граничными условиями.
// Вход: x[0..n-1], f[0..n-1]
// Выход: d[0..n-1] – производные
// ----------------------------------------------------------------------
static void computeDerivatives1D_SignMin(const double* x, const double* f, int n, double* d)
{
    if (n < 2) return;

    // разделённые разности на отрезках
    std::vector<double> dd(n - 1);
    for (int i = 0; i < n - 1; ++i)
        dd[i] = (f[i+1] - f[i]) / (x[i+1] - x[i]);

    if (n == 2) {
        d[0] = d[1] = dd[0];
        return;
    }

    // внутренние узлы (i = 1 .. n-2)
    for (int i = 1; i < n - 1; ++i) {
        double sl = dd[i-1];
        double sr = dd[i];
        if (sl * sr > 0.0) {
            d[i] = (sl > 0.0 ? 1.0 : -1.0) * std::min(std::fabs(sl), std::fabs(sr));
        } else {
            d[i] = 0.0;
        }
    }

    // границы из естественных условий
    // левая: 2*d[0] + d[1] = 3*dd[0]
    d[0] = (3.0 * dd[0] - d[1]) / 2.0;
    // правая: d[n-2] + 2*d[n-1] = 3*dd[n-2]
    d[n-1] = (3.0 * dd[n-2] - d[n-2]) / 2.0;
}

// ----------------------------------------------------------------------
// Кубическая интерполяция на отрезке (та же, что и для 1.13)
// ----------------------------------------------------------------------
static double hermiteInterp1D(double x0, double x1, double y0, double y1, double dy0, double dy1, double x)
{
    double h = x1 - x0;
    double t = (x - x0) / h;
    double t2 = t * t, t3 = t2 * t;

    double H00 = 2*t3 - 3*t2 + 1;
    double H10 = t3 - 2*t2 + t;
    double H01 = -2*t3 + 3*t2;
    double H11 = t3 - t2;

    return y0 * H00 + dy0 * h * H10 + y1 * H01 + dy1 * h * H11;
}

// ----------------------------------------------------------------------
// Вычисление всех производных (экспортируемая)
// ----------------------------------------------------------------------
void computeAllDerivatives31(int nx, int ny,
                             const std::vector<double>& x,
                             const std::vector<double>& y,
                             const std::vector<double>& f_vals,
                             std::vector<double>& dx,
                             std::vector<double>& dy,
                             std::vector<double>& dxy)
{
    dx.assign(nx * ny, 0.0);
    dy.assign(nx * ny, 0.0);
    dxy.assign(nx * ny, 0.0);

    // Производные по x (используем правило с sign и min)
    std::vector<double> f_col(nx);
    std::vector<double> d_col(nx);
    for (int j = 0; j < ny; ++j) {
        
        for (int i = 0; i < nx; ++i)
            f_col[i] = f_vals[i * ny + j];
        
        computeDerivatives1D_SignMin(x.data(), f_col.data(), nx, d_col.data());
        for (int i = 0; i < nx; ++i)
            dx[i * ny + j] = d_col[i];
    }

    // Производные по y
    std::vector<double> f_row(ny);
    std::vector<double> d_row(ny);
    for (int i = 0; i < nx; ++i) {
        
        for (int j = 0; j < ny; ++j)
            f_row[j] = f_vals[i * ny + j];
        
        computeDerivatives1D_SignMin(y.data(), f_row.data(), ny, d_row.data());
        for (int j = 0; j < ny; ++j)
            dy[i * ny + j] = d_row[j];
    }

    // Смешанные производные (∂/∂x от dy)
    std::vector<double> dy_col(nx);
    std::vector<double> dxy_col(nx);
    for (int j = 0; j < ny; ++j) {
        
        for (int i = 0; i < nx; ++i)
            dy_col[i] = dy[i * ny + j];
        
        computeDerivatives1D_SignMin(x.data(), dy_col.data(), nx, dxy_col.data());
        for (int i = 0; i < nx; ++i)
            dxy[i * ny + j] = dxy_col[i];
    }
}

// ----------------------------------------------------------------------
// Вычисление значения интерполянта (экспортируемая)
// ----------------------------------------------------------------------
double GetValue31(double px, double py,
                  const std::vector<double>& x, int nx,
                  const std::vector<double>& y, int ny,
                  const std::vector<double>& f_vals,
                  const std::vector<double>& dx,
                  const std::vector<double>& dy,
                  const std::vector<double>& dxy)
{
    // Поиск отрезка по x
    int ix = std::lower_bound(x.begin(), x.end(), px) - x.begin();
    if (ix == 0) ix = 0;
    else if (ix == nx) ix = nx - 2;
    else ix = ix - 1;
    ix = std::max(0, std::min(nx-2, ix));
    double x0 = x[ix], x1 = x[ix+1];

    // Интерполяция по x для каждого y_j
    std::vector<double> F_j(ny), Fy_j(ny);
    for (int j = 0; j < ny; ++j) {
        double y0 = f_vals[ix * ny + j];
        double y1 = f_vals[(ix+1) * ny + j];
        double dy0 = dx[ix * ny + j];
        double dy1 = dx[(ix+1) * ny + j];
        F_j[j] = hermiteInterp1D(x0, x1, y0, y1, dy0, dy1, px);

        double z0 = dy[ix * ny + j];
        double z1 = dy[(ix+1) * ny + j];
        double dz0 = dxy[ix * ny + j];
        double dz1 = dxy[(ix+1) * ny + j];
        Fy_j[j] = hermiteInterp1D(x0, x1, z0, z1, dz0, dz1, px);
    }

    // Интерполяция по y
    int iy = std::lower_bound(y.begin(), y.end(), py) - y.begin();
    if (iy == 0) iy = 0;
    else if (iy == ny) iy = ny - 2;
    else iy = iy - 1;
    iy = std::max(0, std::min(ny-2, iy));
    double y0 = y[iy], y1 = y[iy+1];
    double v0 = F_j[iy], v1 = F_j[iy+1];
    double dv0 = Fy_j[iy], dv1 = Fy_j[iy+1];

    return hermiteInterp1D(y0, y1, v0, v1, dv0, dv1, py);
}