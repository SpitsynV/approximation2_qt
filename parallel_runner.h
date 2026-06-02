#pragma once
#include <QObject>
#include <QVector>
#include "method_result.h"
#include "shared_input.h"

/// Организатор параллельного запуска всех методов интерполяции.
///
/// Запускает 4 независимых std::thread (по одному на метод),
/// ждёт завершения, применяет MethodSelector и испускает сигнал resultsReady
class ParallelRunner : public QObject {
    Q_OBJECT
public:
    explicit ParallelRunner(QObject* parent = nullptr);

    /// Возвращает заполненный вектор результатов.
    QVector<MethodResult> runAll(const SharedInputData& input, int& bestIdx);

    /// Асинхронный запуск: возвращается немедленно.
    /// По завершении испускает сигнал resultsReady(results, bestIdx).
    void runAllAsync(const SharedInputData& input);

signals:
    /// Испускается из фонового потока; Qt доставляет его в main-thread
    void resultsReady(QVector<MethodResult> results, int bestIdx);
};
