#pragma once
#include <QObject>
#include <QVector>
#include "method_result.h"
#include "shared_input.h"

/// Оркестратор параллельного запуска всех методов интерполяции.
///
/// Запускает 4 независимых std::thread (по одному на метод),
/// ждёт завершения, применяет MethodSelector и испускает сигнал
/// resultsReady с результатами и индексом победителя.
///
/// Гарантии:
/// - SharedInputData копируется внутрь, поэтому Approximator2D
///   можно изменять сразу после вызова runAllAsync().
/// - Сигнал resultsReady всегда приходит в поток, которому принадлежит
///   приёмник (QueuedConnection — стандартное поведение для cross-thread).
class ParallelRunner : public QObject {
    Q_OBJECT
public:
    explicit ParallelRunner(QObject* parent = nullptr);

    /// Синхронный запуск: блокирует текущий поток до завершения всех 4 методов.
    /// Возвращает заполненный (и уже scored) вектор результатов.
    /// bestIdx — индекс победителя (или -1).
    QVector<MethodResult> runAll(const SharedInputData& input, int& bestIdx);

    /// Асинхронный запуск: возвращается немедленно.
    /// По завершении испускает сигнал resultsReady(results, bestIdx).
    void runAllAsync(const SharedInputData& input);

signals:
    /// Испускается из фонового потока; Qt доставляет его в main-thread
    /// через очередь событий (QueuedConnection).
    void resultsReady(QVector<MethodResult> results, int bestIdx);
};
