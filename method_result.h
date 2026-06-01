#pragma once
#include <QString>
#include <QMetaType>
#include <QVector>
/// Результат работы одного метода интерполяции в отдельном потоке.
struct MethodResult {
    int     methodId  = 0;     // 1..4
    QString name;              // название метода
    double  maxError  = 0.0;   // max|approx(x,y) - f(x,y)| на тестовой сетке
    double  elapsedMs = 0.0;   // время выполнения (мс)
    double  score     = 0.0;   // критерий — заполняет MethodSelector
    bool    valid     = false; // true если метод выполнился без ошибок
};

// Регистрируем типы для передачи через очередь сигналов Qt
Q_DECLARE_METATYPE(MethodResult)
Q_DECLARE_METATYPE(QVector<MethodResult>)
