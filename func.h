#pragma once

// f(x,y) по номеру k=0..7 (требование 2)
double GetExactValue(double x, double y, int k);
double GetDX(double x, double y, int k);
double GetDY(double x, double y, int k);
double GetDXY(double x, double y, int k);
// текстовое представление функции
const char *FuncName(int k);