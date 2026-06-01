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