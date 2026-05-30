#include "func.h"
#include <cmath>

double GetExactValue(double x, double y, int k)
{
    switch (k) {
    case 0: return 1.0;
    case 1: return x;
    case 2: return y;
    case 3: return x + y;
    case 4: return sqrt(x * x + y * y);
    case 5: return x * x + y * y;
    case 6: return exp(x * x - y * y);
    case 7: return 1.0 / (25.0 * (x * x + y * y) + 1.0);
    default: return 0.0;
    }
}
double GetDY(double x, double y, int k)
{
    switch (k) {
    case 0: return 0.0;
    case 1: return 0.0;
    case 2: return 1.0;
    case 3: return 1.0;
    case 4: return (y/sqrt(x * x + y * y));
    case 5: return 2.0*y;
    case 6: return (-2.0)*y*exp(x * x - y * y);
    case 7: return (-50.0)*y / ((25.0 * (x * x + y * y) + 1.0)*(25.0 * (x * x + y * y) + 1.0));
    default: return 0.0;
    }
}

double GetDX(double x, double y, int k)
{
    switch (k) {
    case 0: return 0.0;
    case 1: return 1.0;
    case 2: return 0.0;
    case 3: return 1.0;
    case 4: return (x/sqrt(x * x + y * y));
    case 5: return 2.0*x;
    case 6: return (2.0*x*exp(x * x - y * y));
    case 7: return (-50.0*x) / ((25.0 * (x * x + y * y) + 1.0)*(25.0 * (x * x + y * y) + 1.0));
    default: return 0.0;
    }
}
double GetDXY(double x, double y, int k)
{
    switch (k) {
    case 0: return 0.0;
    case 1: return 0.0;
    case 2: return 0.0;
    case 3: return 0.0;
    case 4: return ((-x*y)/pow(x * x + y * y, 1.5));
    case 5: return 0.0;
    case 6: return (-4.0*x*y*exp(x * x - y * y));
    case 7: return (5000.0*x*y) / ((25.0 * (x * x + y * y) + 1.0)*(25.0 * (x * x + y * y) + 1.0)*(25.0 * (x * x + y * y) + 1.0));
    default: return 0.0;
    }
}


const char *FuncName(int k)
{
    switch (k) {
    case 0: return "f(x,y)=1";
    case 1: return "f(x,y)=x";
    case 2: return "f(x,y)=y";
    case 3: return "f(x,y)=x+y";
    case 4: return "f(x,y)=sqrt(x^2+y^2)";
    case 5: return "f(x,y)=x^2+y^2";
    case 6: return "f(x,y)=exp(x^2-y^2)";
    case 7: return "f(x,y)=1/(25*(x^2+y^2)+1)";
    default: return "?";
    }
}