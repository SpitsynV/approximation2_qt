#include "error.h"
#include <vector>
#include <thread>
#include <functional>
#include <cmath>
#include <algorithm> // для std::max


double maxAbsoluteError(double a, double b, double c, double d,
     const std::function<double(double, double)> &exactFunc,
     const std::function<double(double, double)> &approxFunc, int N)
{
    double maxErr = 0.0;
    double stepx = (b - a) / N;
    double stepy= (d-c)/N;
    double diff=0; double x=0; double y=0;
    for (int i = 0; i <= N; i++) {
        x = a + i * stepx;
        for(int j=0;j<=N;j++){
        y=  c+ j*stepy;
        diff = std::abs(exactFunc(x,y) - approxFunc(x,y));
        if (diff > maxErr)
            maxErr = diff;
        }
    }
    return maxErr;
}



double maxAbsoluteErrorParallel(double a, double b, double c, double d,
    const std::function<double(double, double)> &exactFunc,
    const std::function<double(double, double)> &approxFunc,
    int N, int numThreads)
{
    // Если задано 0 потоков или слишком много, используем число hardware потоков
    if (numThreads == 0)
        numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0)
        numThreads = 1;

    // Количество итераций по x
    int totalX = N + 1;
    // Шаги сетки
    double stepx = (b - a) / N;
    double stepy = (d - c) / N;

    // Хранилище для локальных максимумов потоков (размер = numThreads)
    std::vector<double> localMax(numThreads, 0.0);
    std::vector<std::thread> threads;
    threads.reserve(numThreads);

    // Размер базового блока итераций на поток
    int chunkSize = totalX / numThreads;
    int remainder = totalX % numThreads;

    int startIdx = 0;

    for (int t = 0; t < numThreads; ++t) {
        //сколько итераций возьмёт этот поток
        int myCount = chunkSize + (t < remainder ? 1 : 0);
        int endIdx = startIdx + myCount;  // не включая
        //те на блок [StartIdx, endIdx)

        // Запускаем поток
        threads.emplace_back([&, start = startIdx, end = endIdx, t]() {
            double maxErrLocal = 0.0;

            for (int i = start; i < end; ++i) {
                double x = a + i * stepx;
                for (int j = 0; j <= N; ++j) {
                    double y = c + j * stepy;
                    double diff = std::abs(exactFunc(x, y) - approxFunc(x, y));
                    if (diff > maxErrLocal)
                        maxErrLocal = diff;
                }
            }

            // Сохраняем локальный максимум в общий вектор
            localMax[t] = maxErrLocal;
        });

        startIdx = endIdx; // переходим к следующему блоку
    }

    // Ждём завершения всех потоков
    for (auto &th : threads)
        th.join();

    // Финальная редукция – максимум из всех локальных
    double globalMax = 0.0;
    for (double val : localMax)
        if (val > globalMax)
            globalMax = val;

    return globalMax;
}
double integralError2D(
    double a, double b, double c, double d,
    const std::function<double(double,double)> &exactFunc2D,
    const std::function<double(double,double)> &approxFunc2D,
    int Nx, int Ny)  // число узлов по каждой переменной (чётные для Симпсона)
{
    double midX = (a + b) / 2.0, halfX = (b - a) / 2.0;
    double midY = (c + d) / 2.0, halfY = (d - c) / 2.0;
    double hx = M_PI / Nx, hy = M_PI / Ny;
    double sum = 0.0;

    for (int i = 0; i <= Nx; ++i) {
        double theta = i * hx;
        double x = midX + halfX * cos(theta);
        double wx = (i == 0 || i == Nx) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);

        for (int j = 0; j <= Ny; ++j) {
            double phi = j * hy;
            double y = midY + halfY * cos(phi);
            double wy = (j == 0 || j == Ny) ? 1.0 : (j % 2 == 0 ? 2.0 : 4.0);

            double diff = exactFunc2D(x, y) - approxFunc2D(x, y);
            sum += wx * wy * diff * diff;
        }
    }
    return std::sqrt(hx * hy / 9.0 * sum);   // (hx/3)*(hy/3) = hx*hy/9
}