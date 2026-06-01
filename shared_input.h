#pragma once
#include <vector>

//Все данные из approximator2d, передавать в потоки по значению!
struct SharedInputData {
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;  // область [a,b]x[c,d]
    int    nx = 0, ny = 0;    // размер интерполяционной сетки
    int    mx = 0, my = 0;    // размер визуализационной сетки
    int    k  = 0;            // номер функции
    double maxAbsF = 0.0;     // max|f| без возмущения

    std::vector<double> x, y;         //
    std::vector<double> f;            // f[i*ny+j] — значения на сетке (с возмущением)
    std::vector<double> dx, dy, dxy;  //
};
