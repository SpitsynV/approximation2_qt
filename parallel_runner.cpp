#include "parallel_runner.h"
#include "method_runner.h"
#include "method_selector.h"

#include <thread>
#include <vector>

ParallelRunner::ParallelRunner(QObject* parent)
    : QObject(parent)
{}

// ─── Синхронный запуск ────────────────────────────────────────────────────────
QVector<MethodResult> ParallelRunner::runAll(const SharedInputData& input, int& bestIdx)
{
    constexpr int kNumMethods = 4;

    // Создаём runner-ы заранее, чтобы ссылки в лямбдах оставались валидными
    std::vector<MethodRunner> runners;
    runners.reserve(kNumMethods);
    for (int id = 1; id <= kNumMethods; ++id)
        runners.emplace_back(id, input);

    // Запускаем 4 потока
    std::vector<std::thread> threads;
    threads.reserve(kNumMethods);
    for (auto& runner : runners)
        threads.emplace_back([&runner]{ runner.run(); });

    // Ждём завершения всех
    for (auto& t : threads)
        t.join();

    // Собираем результаты
    QVector<MethodResult> results;
    results.reserve(kNumMethods);
    for (const auto& r : runners)
        results.push_back(r.result());

    // Вычисляем score и выбираем победителя
    bestIdx = MethodSelector::selectBest(results);
    return results;
}

// ─── Асинхронный запуск ───────────────────────────────────────────────────────
void ParallelRunner::runAllAsync(const SharedInputData& input)
{
    // Копируем данные — Approximator2D можно менять сразу после возврата
    SharedInputData inputCopy = input;

    // Координирующий поток: запускает 4 рабочих, ждёт, испускает сигнал
    std::thread([this, inputCopy = std::move(inputCopy)]() mutable {
        int bestIdx = -1;
        QVector<MethodResult> results = this->runAll(inputCopy, bestIdx);

        // Вывод в консоль (из любого потока — cout thread-safe)
        MethodSelector::printResults(results, bestIdx);

        // emit из non-main thread: Qt сам доставит через QueuedConnection
        emit resultsReady(results, bestIdx);
    }).detach();
}
