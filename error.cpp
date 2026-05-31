#include "error.h"

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