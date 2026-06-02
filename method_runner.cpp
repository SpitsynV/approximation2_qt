#include "method_runner.h"
#include "method_result.h"

#include "task3.h"    // Бессель
#include "task31.h"   // Разделённые разности
#include "task5.h"    // Эрмит (точные производные)
#include "task6.h"    // МНК Чебышев
#include "func.h"
#include "error.h"

#include <chrono>
#include <cmath>
#include <algorithm>
#include <functional>
#include <stdexcept>

// ─────────────────────────────────────────────────────────────────────────────

MethodRunner::MethodRunner(int methodId, const SharedInputData& input)
    : m_methodId(methodId), m_input(input)
{}

void MethodRunner::run()
{
    switch (m_methodId) {
    case 1: m_result = runMethod1(); break;
    case 2: m_result = runMethod2(); break;
    case 3: m_result = runMethod3(); break;
    case 4: m_result = runMethod4(); break;
    default:
        m_result.valid = false;
        m_result.name  = QStringLiteral("Unknown method %1").arg(m_methodId);
    }
}

// ───расчёт погрешности ────────────────────────────────────────────
//template<typename Fn>
double MethodRunner::computeMaxError(std::function<double(double, double)> approxFn, std::function<double(double, double)> exactFn) const
{
   
    const int tn = std::max(m_input.mx, kMinTestN);
    const int tm = std::max(m_input.my, kMinTestN);
    return maxAbsoluteErrorParallel(m_input.a, m_input.b, m_input.c, m_input.d, exactFn, approxFn, std::max(tn,tm),4);
}

// ─── Метод 1: Бессель (task3) ─────────────────────────────────────────────────
MethodResult MethodRunner::runMethod1()
{
    using Clock = std::chrono::high_resolution_clock;

    MethodResult r;
    r.methodId = 1;
    r.name     = QStringLiteral("Бессель (task3)");
    // Нет гонки-
    // Локальные буферы — не пересекаются с другими потоками
    const int N = m_input.nx * m_input.ny;
    std::vector<double> a3dx(N), a3dy(N), a3dxy(N);

    auto t0 = Clock::now();

    // 1) Построение приближённых производных по методу Бесселя
    computeAllDerivatives(m_input.nx, m_input.ny,
                          m_input.x,  m_input.y, m_input.f,
                          a3dx, a3dy, a3dxy);
    
    // 2) Оценка погрешности на тестовой сетке
    r.maxError = computeMaxError([&](double x, double y) {
        return GetValue13(x, y,
                          m_input.x, m_input.nx,
                          m_input.y, m_input.ny,
                          m_input.f,
                          a3dx, a3dy, a3dxy);
    }, m_input.func);
    auto t1 = Clock::now();
    r.elapsedMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.valid     = true;
    return r;
}

// ─── Метод 2: Разделённые разности (task31) ───────────────────────────────────
MethodResult MethodRunner::runMethod2()
{
    using Clock = std::chrono::high_resolution_clock;

    MethodResult r;
    r.methodId = 2;
    r.name     = QStringLiteral("Разд. разности (task31)");

    const int N = m_input.nx * m_input.ny;
    std::vector<double> a4dx(N), a4dy(N), a4dxy(N);

    auto t0 = Clock::now();

    computeAllDerivatives31(m_input.nx, m_input.ny,
                            m_input.x,  m_input.y, m_input.f,
                            a4dx, a4dy, a4dxy);
    
    r.maxError = computeMaxError([&](double x, double y) {
        return GetValue31(x, y,
                          m_input.x, m_input.nx,
                          m_input.y, m_input.ny,
                          m_input.f,
                          a4dx, a4dy, a4dxy);
    }, m_input.func);

    auto t1 = Clock::now();
    r.elapsedMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.valid     = true;
    return r;
}

// ─── Метод 3: Эрмит, точные производные (task5) ───────────────────────────────
MethodResult MethodRunner::runMethod3()
{
    using Clock = std::chrono::high_resolution_clock;

    MethodResult r;
    r.methodId = 3;
    r.name     = QStringLiteral("Эрмит (task5)");

    // Точные производные уже вычислены в SharedInputData::dx/dy/dxy —
    // этот метод не требует этапа построения коэффициентов.

    auto t0 = Clock::now();

    r.maxError =computeMaxError([&](double x, double y) {
        return GetValue5(x, y,
                          m_input.x, m_input.nx,
                          m_input.y, m_input.ny,
                          m_input.f,
                          m_input.dx, m_input.dy, m_input.dxy);
    }, m_input.func);

    auto t1 = Clock::now();
    r.elapsedMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.valid     = true;
    return r;
}

// ─── Метод 4: МНК по Чебышеву (task6) ────────────────────────────────────────
MethodResult MethodRunner::runMethod4()
{
    using Clock = std::chrono::high_resolution_clock;

    MethodResult r;
    r.methodId = 4;
    r.name     = QStringLiteral("МНК Чебышев (task6)");

    std::vector<double> c6(m_input.nx * m_input.ny, 0.0);

    auto t0 = Clock::now();

    // Построение коэффициентов: GetCoefficients6 принимает точную функцию
    
    GetCoefficients6(m_input.nx, m_input.ny,
                     m_input.a,  m_input.b,
                     m_input.c,  m_input.d,
                     m_input.func, c6);

    r.maxError = computeMaxError([&](double x, double y) {
        return GetValue6(x, y,
                         m_input.a, m_input.b,
                         m_input.c, m_input.d,
                         c6,
                         m_input.nx, m_input.ny);
    }, m_input.func);

    auto t1 = Clock::now();
    r.elapsedMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    r.valid     = true;
    return r;
}
