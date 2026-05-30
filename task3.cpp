// task3.cpp
// Реализация задачи 1.13: кусочная интерполяция кубическими многочленами Эрмита
// с естественными граничными условиями. Алгоритм основан на последовательной
// одномерной интерполяции (тензорное произведение), без явного построения
// бикубических коэффициентов. Соответствует Главе II, §1.4 и §1.6.3.
#include"task3.h"
#include <vector>
#include <cmath>
#include <algorithm>

// ----------------------------------------------------------------------
// Одномерное вычисление производных d[i] для сетки x[0..n-1] и значений f[0..n-1]
// с использованием естественных граничных условий (P''=0 на концах).
// Формулы: (7) из §19.4 (стр. 69) и условия из §19.8.3 (стр. 76)
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

    // внутренние узлы: среднее арифметическое (метод Бесселя)
    for (int i = 1; i < n - 1; ++i)
        d[i] = 0.5 * (dd[i-1] + dd[i]);

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

    // базисные функции Эрмита (стр. 68 неявно, но можно взять из общей теории)
    double H00 = 2*t3 - 3*t2 + 1;   // вклад y0
    double H10 = t3 - 2*t2 + t;     // вклад dy0 * h
    double H01 = -2*t3 + 3*t2;      // вклад y1
    double H11 = t3 - t2;           // вклад dy1 * h

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
    dx.resize(nx * ny);
    dy.resize(nx * ny);
    dxy.resize(nx * ny);

    // 1. Производные по x (dx) для каждой строки j
    for (int j = 0; j < ny; ++j) {
        std::vector<double> f_col(nx);
        for (int i = 0; i < nx; ++i)
            f_col[i] = f_vals[i * ny + j];
        std::vector<double> d_col(nx);
        computeDerivatives1D(x.data(), f_col.data(), nx, d_col.data());
        for (int i = 0; i < nx; ++i)
            dx[i * ny + j] = d_col[i];
    }

    // 2. Производные по y (dy) для каждого столбца i
    for (int i = 0; i < nx; ++i) {
        std::vector<double> f_row(ny);
        for (int j = 0; j < ny; ++j)
            f_row[j] = f_vals[i * ny + j];
        std::vector<double> d_row(ny);
        computeDerivatives1D(y.data(), f_row.data(), ny, d_row.data());
        for (int j = 0; j < ny; ++j)
            dy[i * ny + j] = d_row[j];
    }

    // 3. Смешанные производные dxy: применяем computeDerivatives1D к dy по x
    for (int j = 0; j < ny; ++j) {
        std::vector<double> dy_col(nx);
        for (int i = 0; i < nx; ++i)
            dy_col[i] = dy[i * ny + j];
        std::vector<double> dxy_col(nx);
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

// ----------------------------------------------------------------------
// Основная функция построения коэффициентов (для совместимости с интерфейсом).
// В данном алгоритме коэффициенты не хранятся в явном виде, поэтому эта функция
// только вычисляет производные и сохраняет их в глобальные массивы (или возвращает).
// Однако для удобства использования можно вычислить производные и сохранить их,
// а для вычисления значения использовать GetValue13 с этими производными.
// Вместо coeffs мы возвращаем производные через дополнительные векторы.
// В вашем проекте вы можете адаптировать вызов: сначала вызвать computeAllDerivatives,
// затем при вычислении значения использовать GetValue13, передавая туда производные.
// ----------------------------------------------------------------------
