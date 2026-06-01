#include "task6.h"
#include <cmath>
#include <vector>
#include <algorithm>

// ----------------------------------------------------------------------
// Одномерное разложение по многочленам Чебышева (метод наименьших квадратов).
// Используется то что было в aprox1_qt
// Вход:
//   n   – число узлов
//   a,b – отрезок [a,b]
//   f   – значения функции в узлах Чебышева (размер n)
// Выход:
//   coef – коэффициенты α_i (размер n)
// ----------------------------------------------------------------------
static void chebyshev1D(int n, double a, double b,
                        const std::vector<double>& f,
                        std::vector<double>& coef)
{
    // Генерируем узлы Чебышева на [a,b] в порядке возрастания.
    std::vector<double> x(n);
    for (int k = 0; k < n; ++k) {
        double t = cos(M_PI * (2.0 * (n - 1 - k) + 1.0) / (2.0 * n));
        x[k] = (a + b) / 2.0 + (b - a) / 2.0 * t;
    }

    coef.assign(n, 0.0);
    double gij = 0.0, val = 0.0;
    double zj = 0.0;
    std::vector<double> gk0(n, 0.0);
    std::vector<double> gk1(n, 0.0);
    std::vector<double> z(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == 0) {
                gij = f[j];
                gk0[j] = gij;
            } else if (i == 1) {
                double t = (2.0 * x[j] - (a + b)) / (b - a); // t ∈ [-1,1]
                zj = t;
                gij = zj * f[j];
                gk1[j] = gij;
                z[j] = 2.0 * zj;
            } else {
                if (i % 2 == 0) {
                    gij = z[j] * gk1[j] - gk0[j];
                    gk0[j] = gij;
                } else {
                    gij = z[j] * gk0[j] - gk1[j];
                    gk1[j] = gij;
                }
            }
            val += gij;
        }
        if (i == 0)
            coef[i] = val / n;
        else
            coef[i] = 2.0 * val / n;
        val = 0.0;
    }
}

// ----------------------------------------------------------------------
// Двумерное тензорное произведение (использует одномерный метод)
// ----------------------------------------------------------------------
void GetCoefficients6(int nx, int ny,
                      double a, double b, double c, double d,
                      std::function<double(double, double)> func,
                      std::vector<double>& coeffs)
{
    // 1. Генерируем узлы Чебышева по x и y (в том же порядке, что в chebyshev1D)
    std::vector<double> x_cheb(nx), y_cheb(ny);
    double t=0;
    for (int i = 0; i < nx; ++i) {
        t  = cos(M_PI * (2.0*(nx-1-i) + 1.0) / (2.0*nx));
        x_cheb[i] = (a + b)/2 + (b - a)/2 * t;
    }
    for (int j = 0; j < ny; ++j) {
        t = cos(M_PI * (2.0 * (ny - 1 - j) + 1.0) / (2.0 * ny));
        y_cheb[j] = (c + d)/2 + (d - c)/2 * t;
    }
    
    // 2. Вычисляем значения функции в узлах сетки
    std::vector<double> f_vals(nx * ny);
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            f_vals[i * ny + j] = func(x_cheb[i], y_cheb[j]);
        }
    }

    // 3. Промежуточная матрица temp (размер nx*ny)
    std::vector<double> temp(nx * ny);
    std::vector<double> col(nx), alpha_col(nx);
    for (int j = 0; j < ny; ++j) {
        // Забираем столбец j (фиксированный y) – значения f(x_i, y_j)
        for (int i = 0; i < nx; ++i)
            col[i] = f_vals[i * ny + j];
        // Одномерное преобразование по x
        chebyshev1D(nx, a, b, col, alpha_col);
        for (int i = 0; i < nx; ++i)
            temp[i * ny + j] = alpha_col[i];
    }

    // 4. Итоговые коэффициенты (размер nx*ny)
    coeffs.resize(nx * ny);
    std::vector<double> row(ny), alpha_row(ny);
    for (int i = 0; i < nx; ++i) {
        // Забираем строку i (фиксированный x) из temp
        for (int j = 0; j < ny; ++j)
            row[j] = temp[i * ny + j];
        // Одномерное преобразование по y
        chebyshev1D(ny, c, d, row, alpha_row);
        for (int j = 0; j < ny; ++j)
            coeffs[i * ny + j] = alpha_row[j];
    }
}

// ----------------------------------------------------------------------
// Вычисление значения разложения в точке (x,y)
// ----------------------------------------------------------------------
double GetValue6(double x, double y,
                 double a, double b, double c, double d,
                 const std::vector<double>& coeffs,
                 int nx, int ny)
{
    // Приведение координат к отрезку [-1,1]
    double xi = 2.0 * (x - a) / (b - a) - 1.0;
    double eta = 2.0 * (y - c) / (d - c) - 1.0;
    xi = std::max(-1.0, std::min(1.0, xi));
    eta = std::max(-1.0, std::min(1.0, eta));

    // Предвычисление многочленов Чебышева T_i(xi) и T_j(eta)
    std::vector<double> Txi(nx), Teta(ny);
    auto cheb = [](int n, double t) -> double {
        if (n == 0) return 1.0;
        if (n == 1) return t;
        double t0 = 1.0, t1 = t;
        for (int k = 2; k <= n; ++k) {
            double t2 = 2.0 * t * t1 - t0;
            t0 = t1;
            t1 = t2;
        }
        return t1;
    };
    for (int i = 0; i < nx; ++i) Txi[i] = cheb(i, xi);
    for (int j = 0; j < ny; ++j) Teta[j] = cheb(j, eta);

    double result = 0.0;
    for (int i = 0; i < nx; ++i) {
        double row_sum = 0.0;
        const double* row = &coeffs[i * ny];
        for (int j = 0; j < ny; ++j)
            row_sum += row[j] * Teta[j];
        result += Txi[i] * row_sum;
    }
    return result;
}