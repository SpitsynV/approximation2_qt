#pragma once
#include "method_result.h"
#include <QVector>
#include <QString>

/// Метод с наименьшим score выигрывает.
class MethodSelector {
public:
    static constexpr double kAlpha = 0.7;
    static int selectBest(QVector<MethodResult>& results);
    static void printResults(const QVector<MethodResult>& results, int bestIdx);
    static QString toHtmlTable(const QVector<MethodResult>& results, int bestIdx);
};
