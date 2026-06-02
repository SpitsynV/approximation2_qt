#pragma once
#include <vector>
#include <functional>
#include <QString>
#include <QColor>
#include "error.h"
#include "shared_input.h"

// Аналог Approximator из 1D-задачи — управляет состоянием и пересчётом
class Approximator2D
{
public:
    Approximator2D(double a, double b, double c, double d,
                   int nx, int ny, int mx, int my, int k);

    // Пересчёт коэффициентов при изменении параметров
    void rebuild();

    // Значение исходной функции (с возмущением) и всех приближений
    double f   (double x, double y) const;
    double approx1(double x, double y) const;   //Бессель
    double approx2(double x, double y) const;   //Раздел разности
    double approx3(double x, double y) const;   //Эрмит
    double approx4(double x, double y) const;   //МНК Чебышев

    double getMaxError1() const;
    double getMaxError2() const;
    double getMaxError3() const;
    double getError4() const;

    // Геттеры
    double a() const { return m_a; }
    double b() const { return m_b; }
    double c() const { return m_c; }
    double d() const { return m_d; }
    int nx() const { return m_nx; }
    int ny() const { return m_ny; }
    int mx() const { return m_mx; }
    int my() const { return m_my; }
    int k()  const { return m_k;  }
    int p()  const { return m_p;  }
    int scaleExp()  const { return m_scale; }
    int graphMode() const { return m_graphMode; }
    double angle()  const { return m_angle; }   // угол поворота вокруг OZ (рад)

    // Сеттеры
    void setNx(int nx);
    void setNy(int ny);
    void setP(int p);
    void setScaleExp(int s) { m_scale = s; }
    void setGraphMode(int m){ m_graphMode = m; }
    void setAngle(double a) { m_angle = a; }

    // Циклический перебор
    void nextK();
    void nextGraphMode();

    QString functionName() const;

    // Узлы интерполяции
    const std::vector<double>& getX()     const { return m_x; }
    const std::vector<double>& getY()     const { return m_y; }
    const std::vector<double>& getFvals() const { return m_f; }

    // Возвращает функцию для отрисовки (с учётом graphMode)
    // functor(x,y) -> значение
    std::function<double(double, double)> getPlotFunc() const;
    QString getPlotName() const;
    //
    SharedInputData extractSharedInput() const;
    ///Результат для лучшего метода
    int IDX=-1;
    void setIDX(int i){IDX=i;};




private:
    double m_a, m_b, m_c, m_d;   // область [a,b]x[c,d]
    int    m_nx, m_ny;            // число точек интерполяции
    int    m_mx, m_my;            // число точек визуализации
    int    m_k;                   // номер функции
    int    m_p;                   // возмущение
    int    m_scale;               // показатель масштаба s
    int    m_graphMode;           // 0..4
    double m_angle;               // угол поворота вокруг OZ

    double m_maxAbsF;             // max|f| без возмущения

    //сетка интерполяции
    std::vector<double> m_x, m_y;
    std::vector<double> m_f;      // m_f[i*m_ny+j] = f(m_x[i], m_y[j]) с возмущением
    std::vector<double> m_dx;     // df/dx по той же формуле(точные)
    std::vector<double> m_dy;
    std::vector<double> m_dxy;    // d^2f/dxdy(точные)
    // Коэффициенты метода 3
    //std::vector<double> m_c3;
    std::vector<double> m_a3dx;     // для метода 3 апрокс производные
    std::vector<double> m_a3dy;
    std::vector<double> m_a3dxy; 

    // Коэффициенты метода 31 (кусочный кубический)
    //std::vector<double> m_c4;
    std::vector<double> m_a4dx;     // для метода 4 апрокс производные
    std::vector<double> m_a4dy;
    std::vector<double> m_a4dxy;
    // Коэффициенты метода 5
    //std::vector<double> m_c5;
    /*Не нужны производные-Уже посчитаны точные производные*/

    //Коэф метода 6
    std::vector<double> m_c6; //Коэф метода 6
                            //Метод не требуе производных

    void initGrid();
    
};