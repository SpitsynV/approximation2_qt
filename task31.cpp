#include "task31.h"
#include <algorithm>
#include <cmath>

// ─────────────────────────────────────────────────────────
// 1D-вспомогательные функции (§19.1, §19.6, §19.8.3)
// ─────────────────────────────────────────────────────────

// Разделённая разность первого порядка f(x_i; x_{i+1}) (§I.4, ф.(4.1))
static inline double DD1(double fi, double fi1, double xi, double xi1)
{
    return (fi1 - fi) / (xi1 - xi);
}

// Вычисление параметров d[0..n-1] по §19.6 (ф.(9) стр.71) + §19.8.3 (стр.76)
// vals[i] = значение функции в x[i],  nodes[i] = x[i]
// Индексация 0-based; в книге i = 1..n
static void ComputeD_DD_NaturalBC(int n, const double *nodes,
                                   const double *vals, double *d)
{
    if (n == 1) { d[0] = 0.0; return; }
    if (n == 2) {
        // только два узла — линейная интерполяция
        double r = DD1(vals[0], vals[1], nodes[0], nodes[1]);
        d[0] = r;
        d[1] = r;
        return;
    }

    // Вычисляем разделённые разности первого порядка
    std::vector<double> dd(n - 1);
    for (int i = 0; i < n - 1; i++)
        dd[i] = DD1(vals[i], vals[i + 1], nodes[i], nodes[i + 1]);

    // §19.6, ф.(9): внутренние узлы i=1..n-2 (0-based)
    for (int i = 1; i <= n - 2; i++) {
        double s0 = std::copysign(1.0, dd[i - 1]);
        double s1 = std::copysign(1.0, dd[i]);
        if (s0 == s1)
            d[i] = s1 * std::min(std::fabs(dd[i - 1]), std::fabs(dd[i]));
        else
            d[i] = 0.0;
    }

    // §19.8.3: естественные ГУ — (Pf)''(x_1) = 0, (Pf)''(x_n) = 0
    //   2*d_1 + d_2   = 3*f(x_1;x_2)       (книга: индексы 1-based)
    //   d_{n-1}+2*d_n = 3*f(x_{n-1};x_n)
    d[0]     = (3.0 * dd[0]     - d[1])     / 2.0;
    d[n - 1] = (3.0 * dd[n - 2] - d[n - 2]) / 2.0;
}

// Коэффициенты кубического многочлена P_i на [x_i, x_{i+1}] (§19.2, ф.(3) стр.68):
//   P_i(x) = c1 + c2*(x-xi) + c3*(x-xi)^2 + c4*(x-xi)^3
//   c1 = f(x_i)
//   c2 = d_i
//   c3 = (3*f(x_i;x_{i+1}) - 2*d_i - d_{i+1}) / (x_{i+1}-x_i)
//   c4 = (d_i + d_{i+1} - 2*f(x_i;x_{i+1})) / (x_{i+1}-x_i)^2
static void CubicCoef(double fi, double fi1, double di, double di1,
                      double h, double c[4])
{
    double r = DD1(fi, fi1, 0.0, h); // f(x_i;x_{i+1}) при h = x_{i+1}-x_i
    c[0] = fi;
    c[1] = di;
    c[2] = (3.0 * r - 2.0 * di - di1) / h;
    c[3] = (di + di1 - 2.0 * r) / (h * h);
}

// ─────────────────────────────────────────────────────────
// Матрица A^{-1}(h) (§1.6.3 стр.106)
// Преобразует вектор (f(x_i), d_i, f(x_{i+1}), d_{i+1})^T
// в вектор коэффициентов (c1,c2,c3,c4)^T — то есть CubicCoef
//
// Отдельная функция умножения Gamma = A^{-1}(hi)*F*(A^{-1}(hj))^T
// Использует прямое применение CubicCoef для каждого столбца/строки
// ─────────────────────────────────────────────────────────

// Применяем A^{-1}(h) к каждому столбцу матрицы 4×4 M (in-place по строкам)
// т.е. заменяем M -> A^{-1}(h) * M
// Входной M[r][c]: r — строка (0..3), c — столбец (0..3)
// После: M[r][c] содержит r-ю компоненту коэффициентов для c-го столбца
static void ApplyAinvLeft(double M[4][4], double h)
{
    // каждый столбец c: вход = [M[0][c], M[1][c], M[2][c], M[3][c]]
    //                   = [f_i, d_i, f_{i+1}, d_{i+1}]
    // выход = [c1, c2, c3, c4] из CubicCoef
    for (int c = 0; c < 4; c++) {
        double coef[4];
        CubicCoef(M[0][c], M[2][c], M[1][c], M[3][c], h, coef);
        M[0][c] = coef[0];
        M[1][c] = coef[1];
        M[2][c] = coef[2];
        M[3][c] = coef[3];
    }
}

// Применяем A^{-1}(h) к каждой строке матрицы 4×4 M (умножение справа на (A^{-1})^T)
// т.е. заменяем M -> M * (A^{-1}(h))^T
static void ApplyAinvRight(double M[4][4], double h)
{
    // каждая строка r: вход = [M[r][0], M[r][1], M[r][2], M[r][3]]
    //                        = [f_j, e_j, f_{j+1}, e_{j+1}]
    for (int r = 0; r < 4; r++) {
        double coef[4];
        CubicCoef(M[r][0], M[r][2], M[r][1], M[r][3], h, coef);
        M[r][0] = coef[0];
        M[r][1] = coef[1];
        M[r][2] = coef[2];
        M[r][3] = coef[3];
    }
}

// ─────────────────────────────────────────────────────────
// Основные функции метода 31
// ─────────────────────────────────────────────────────────

void GetCoefficients31(int nx, int ny,
                       const std::vector<double> &x,
                       const std::vector<double> &y,
                       const std::vector<double> &f_vals,
                       std::vector<double> &coef)
{
    // f[i*ny+j] = f(x[i], y[j])

    // ── Шаг 1: x-производные dx[i*ny+j] = d_i^x(y_j) ──────────────
    // Для каждого фиксированного j применяем §19.6+§19.8.3 по x
    std::vector<double> dx(nx * ny), dy(nx * ny), dxy(nx * ny);

    {
        std::vector<double> col(nx);
        std::vector<double> dcol(nx);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) col[i] = f_vals[i * ny + j];
            ComputeD_DD_NaturalBC(nx, x.data(), col.data(), dcol.data());
            for (int i = 0; i < nx; i++) dx[i * ny + j] = dcol[i];
        }
    }

    // ── Шаг 2: y-производные dy[i*ny+j] = e_j^y(x_i) ──────────────
    // Для каждого фиксированного i применяем §19.6+§19.8.3 по y
    {
        std::vector<double> row(ny);
        std::vector<double> drow(ny);
        for (int i = 0; i < nx; i++) {
            for (int j = 0; j < ny; j++) row[j] = f_vals[i * ny + j];
            ComputeD_DD_NaturalBC(ny, y.data(), row.data(), drow.data());
            for (int j = 0; j < ny; j++) dy[i * ny + j] = drow[j];
        }
    }

    // ── Шаг 3: смешанные производные dxy[i*ny+j] = d_i^x(e_j^y(f)) ─
    // Применяем §19.6+§19.8.3 по x к массиву dy[·][j] для каждого j
    {
        std::vector<double> col(nx);
        std::vector<double> dcol(nx);
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) col[i] = dy[i * ny + j];
            ComputeD_DD_NaturalBC(nx, x.data(), col.data(), dcol.data());
            for (int i = 0; i < nx; i++) dxy[i * ny + j] = dcol[i];
        }
    }

    // ── Шаг 4: коэффициенты ячеек Gamma = A^{-1}(h_i)*F*( A^{-1}(h_j) )^T ─
    // (§1.6.3, стр.105-106)
    const int ncx = nx - 1;  // число интервалов по x
    const int ncy = ny - 1;  // число интервалов по y
    coef.resize(ncx * ncy * 16);

    for (int i = 0; i < ncx; i++) {
        double hi = x[i + 1] - x[i];

        for (int j = 0; j < ncy; j++) {
            double hj = y[j + 1] - y[j];

            // Матрица F_{ij} (ф. стр.106):
            //   строки соответствуют k=1..4: (f_i, d_i^x, f_{i+1}, d_{i+1}^x) применённым к y-функционалу
            //   столбцы соответствуют l=1..4: (f_j, e_j^y, f_{j+1}, e_{j+1}^y)
            double F[4][4];
            F[0][0] = f_vals[i       * ny + j];
            F[0][1] = dy    [i       * ny + j];
            F[0][2] = f_vals[i       * ny + j + 1];
            F[0][3] = dy    [i       * ny + j + 1];

            F[1][0] = dx    [i       * ny + j];
            F[1][1] = dxy   [i       * ny + j];
            F[1][2] = dx    [i       * ny + j + 1];
            F[1][3] = dxy   [i       * ny + j + 1];

            F[2][0] = f_vals[(i + 1) * ny + j];
            F[2][1] = dy    [(i + 1) * ny + j];
            F[2][2] = f_vals[(i + 1) * ny + j + 1];
            F[2][3] = dy    [(i + 1) * ny + j + 1];

            F[3][0] = dx    [(i + 1) * ny + j];
            F[3][1] = dxy   [(i + 1) * ny + j];
            F[3][2] = dx    [(i + 1) * ny + j + 1];
            F[3][3] = dxy   [(i + 1) * ny + j + 1];

            // Gamma = A^{-1}(h_i) * F * (A^{-1}(h_j))^T
            ApplyAinvLeft (F, hi);   // F <- A^{-1}(hi) * F
            ApplyAinvRight(F, hj);   // F <- F * (A^{-1}(hj))^T

            // сохраняем gamma[k][l] = F[k][l]
            int base = (i * ncy + j) * 16;
            for (int k = 0; k < 4; k++)
                for (int l = 0; l < 4; l++)
                    coef[base + k * 4 + l] = F[k][l];
        }
    }
}

// Вычисление значения в точке (px, py):
// 1. Бинарный поиск ячейки i (x[i] <= px < x[i+1]), j (y[j] <= py < y[j+1])
// 2. Вычисление кубического двумерного многочлена вложенной схемой Горнера:
//    P(x,y) = sum_{k=0..3} sum_{l=0..3} gamma[k][l] * dx^k * dy^l
//    внутренний Горнер по l (для каждого k): p_k(dy)
//    внешний Горнер по k:  sum_k p_k(dy) * dx^k
double GetValue31(double px, double py,
                  const std::vector<double> &x, int nx,
                  const std::vector<double> &y, int ny,
                  const std::vector<double> &coef)
{
    const int ncx = nx - 1;
    const int ncy = ny - 1;

    // бинарный поиск по x
    int i = (int)(std::upper_bound(x.begin(), x.end(), px) - x.begin()) - 1;
    if (i < 0)   i = 0;
    if (i >= ncx) i = ncx - 1;

    // бинарный поиск по y
    int j = (int)(std::upper_bound(y.begin(), y.end(), py) - y.begin()) - 1;
    if (j < 0)   j = 0;
    if (j >= ncy) j = ncy - 1;

    double dx = px - x[i];
    double dy = py - y[j];

    const double *g = &coef[(i * ncy + j) * 16];  // g[k*4+l] = gamma_{k,l}

    // вложенная схема Горнера: P = sum_k (sum_l g[k][l]*dy^l) * dx^k
    double result = 0.0;
    for (int k = 3; k >= 0; k--) {
        double pk = 0.0;
        for (int l = 3; l >= 0; l--) {
            pk = pk * dy + g[k * 4 + l];
        }
        result = result * dx + pk;
    }
    return result;
}