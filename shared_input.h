#pragma once
#include <vector>
#include<functional>
//Все данные из approximator2d, передавать в потоки по значению!
struct SharedInputData {
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    int    nx = 0, ny = 0;    
    int    mx = 0, my = 0;    
    int    k  = 0;            
    double maxAbsF = 0.0;
    std::function<double(double, double)> func;   // f(x,y) с возмущением
    std::vector<double> x, y;         //
    std::vector<double> f;            // f[i*ny+j] — значения
    std::vector<double> dx, dy, dxy;  //
};
