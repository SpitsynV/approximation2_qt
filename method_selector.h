#pragma once
#include "method_result.h"
#include <QVector>
#include <QString>

/// Вычисляет интегральный критерий качества для каждого метода
/// и определяет наилучший.
///
/// Критерий (меньше = лучше):
///   score_i = α · (E_i / E_max) + (1-α) · (T_i / T_max)
///
/// где E_i — максимальная абсолютная погрешность метода i,
///     T_i — время выполнения метода i (мс),
///     α   — вес погрешности (kAlpha = 0.7 по умолчанию).
///
/// Смысл: нормируем оба показателя в [0,1] и берём взвешенную сумму.
/// Метод с наименьшим score выигрывает.
class MethodSelector {
public:
    /// Вес погрешности в критерии (0 < kAlpha < 1).
    /// Увеличить — приоритет точности; уменьшить — приоритет скорости.
    static constexpr double kAlpha = 0.7;

    /// Заполняет поле MethodResult::score для каждого результата.
    /// Возвращает индекс лучшего элемента в векторе results.
    static int selectBest(QVector<MethodResult>& results);

    /// Выводит таблицу результатов + победителя в stdout (std::cout).
    static void printResults(const QVector<MethodResult>& results, int bestIdx);

    /// Возвращает HTML-таблицу для отображения в QLabel (setTextFormat(Qt::RichText)).
    static QString toHtmlTable(const QVector<MethodResult>& results, int bestIdx);
};
