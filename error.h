#ifndef ERROR_H
#define ERROR_H

#include <vector>
#include <functional>
#include<cmath>
#include <functional>

double maxAbsoluteError(double a, double b, double c, double d,
     const std::function<double(double, double)> &exactFunc,
     const std::function<double(double, double)> &approxFunc, int N=1000);
     
double maxAbsoluteErrorParallel(double a, double b, double c, double d,
    const std::function<double(double, double)> &exactFunc,
    const std::function<double(double, double)> &approxFunc,
    int N=1000, int numThreads=4);

double integralError2D(
    double a, double b, double c, double d,
    const std::function<double(double,double)> &exactFunc2D,
    const std::function<double(double,double)> &approxFunc2D,
    int Nx = 100, int Ny = 100);
#endif