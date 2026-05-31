#include "approximator2d.h"
#include "func.h"
#include "task3.h"
#include "task31.h"
#include <cmath>
#include <algorithm>

Approximator2D::Approximator2D(double a, double b, double c, double d,
                               int nx, int ny, int mx, int my, int k)
    : m_a(a), m_b(b), m_c(c), m_d(d)
    , m_nx(nx), m_ny(ny), m_mx(mx), m_my(my)
    , m_k(k), m_p(0), m_scale(0), m_graphMode(0), m_angle(0.0)
    , m_maxAbsF(0.0)
{
    m_x.resize(m_nx);
    m_y.resize(m_ny);
    m_f.resize(m_nx*m_ny);
    m_dx.resize(m_nx*m_ny);
    m_dy.resize(m_nx*m_ny);
    m_dxy.resize(m_nx*m_ny);
    m_adx.resize(m_nx*m_ny);
    m_ady.resize(m_nx*m_ny);
    m_adxy.resize(m_nx*m_ny);
    m_c3.resize(16*(m_nx-1)*(m_ny-1), 0.0);
    m_c31.resize(16*(m_nx-1)*(m_ny-1),0.0);
    initGrid();
}

void Approximator2D::initGrid()
{
    // сетка интерполяции по x и y.
    m_x.resize(m_nx);
    m_y.resize(m_ny);
    for (int i = 0; i < m_nx; i++){
        m_x[i] = m_a + i * (m_b - m_a) / (m_nx - 1);
    }
    for (int j = 0; j < m_ny; j++){
        m_y[j] = m_c + j * (m_d - m_c) / (m_ny - 1);
    }
    // значения функции//
    m_f.resize(m_nx * m_ny);
    m_maxAbsF = 0.0;
    for (int i = 0; i < m_nx; i++) {
        for (int j = 0; j < m_ny; j++) {
            double v = GetExactValue(m_x[i], m_y[j], m_k);
            m_f[i * m_ny + j] = v;
            //m_dx[i*m_ny+j]= GetDX(m_x[i], m_y[j], m_k);
            //m_dy[i*m_ny+j]= GetDY(m_x[i], m_y[j], m_k);
            //m_dxy[i*m_ny+j]= GetDXY(m_x[i], m_y[j], m_k);
            m_maxAbsF = std::max(m_maxAbsF, std::fabs(v));
        }
    }

    // возмущение: добавляем p*0.1*max|f| к f(x_{nx/2}, y_{ny/2})
    if (m_p != 0) {
        int imid = m_nx / 2;
        int jmid = m_ny / 2;
        m_f[imid * m_ny + jmid] += m_p * 0.1 * m_maxAbsF;
    }
}

void Approximator2D::rebuild()
{
    initGrid();
    computeAllDerivatives(m_nx,m_ny, m_x,m_y,m_f,m_adx, m_ady,m_adxy);
    computeAllDerivatives31(m_nx,m_ny, m_x,m_y,m_f,m_dx, m_dy,m_dxy);
}

void Approximator2D::setNx(int nx) { m_nx = std::max(2, nx); }
void Approximator2D::setNy(int ny) { m_ny = std::max(2, ny); }

void Approximator2D::nextK()
{
    m_k = (m_k + 1) % 8;  // k=0..7
}

void Approximator2D::nextGraphMode()
{
    m_graphMode = (m_graphMode + 1) % 5;  // 0..4
}

QString Approximator2D::functionName() const
{
    return QString(FuncName(m_k));
}

// исходная функция с возмущением
double Approximator2D::f(double x, double y) const
{
    double v = GetExactValue(x, y, m_k);
    // возмущение применяется только в 1 точке
    if (m_p != 0 && m_nx > 0 && m_ny>0) {
        int midx = m_nx / 2;
        int midy = m_ny/2;
        // Сравнение с mid с учётом погрешности
        if (std::abs(x - m_x[midx]) < 1e-12 * std::max(1.0, std::abs(x))) {
            if(std::abs(y - m_y[midy]) < 1e-12 * std::max(1.0, std::abs(y))){
                v += m_p * 0.1 * m_maxAbsF;
            }
        }
    }
    return v;
}

double Approximator2D::approx1(double x, double y) const
{
    return GetValue13(x,y,m_x,m_nx,m_y, m_ny,m_f,m_adx,m_ady,m_adxy);
}

double Approximator2D::approx2(double x, double y) const
{
    return GetValue31(x, y,m_x,m_nx,m_y, m_ny,m_f,m_dx,m_dy,m_dxy );
}
std::function<double(double, double)> Approximator2D::getPlotFunc() const
{
    switch (m_graphMode) {
    case 0:  // функция
        return [this](double x, double y){ return this->f(x, y); };
    case 1:  // приближение по методу 1 
        return [this](double x, double y){ return this->approx1(x, y); };
    case 2:  // погрешность метода 1
        return [this](double x, double y){ return this->approx1(x, y) - this->f(x, y); };
    case 3:  // приближение по методу 2 
        return [this](double x, double y){ return this->approx2(x, y); };
    case 4:  // погрешность метода 2
        return [this](double x, double y){ return this->approx2(x, y) - this->f(x, y); };
    default:
        return [](double, double){ return 0.0; };
    }
}

QString Approximator2D::getPlotName() const
{
    switch (m_graphMode) {
    case 0: return "f(x,y)";
    case 1: return "Approx1 (Newton)";
    case 2: return "Error1 = P1-f";
    case 3: return "Approx2 (PieceCubic)";
    case 4: return "Error2 = P2-f";
    default: return "";
    }
}