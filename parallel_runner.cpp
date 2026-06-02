#include "parallel_runner.h"
#include "method_runner.h"
#include "method_selector.h"

#include <thread>
#include <vector>

ParallelRunner::ParallelRunner(QObject* parent)
    : QObject(parent)
{}

QVector<MethodResult> ParallelRunner::runAll(const SharedInputData& input, int& bestIdx)
{
    const int kNumMethods = 4;

    //Только так
    std::vector<MethodRunner> runners;
    runners.reserve(kNumMethods);
    for (int id = 1; id <= kNumMethods; ++id)
        runners.emplace_back(id, input);

    // Запуск потоков
    std::vector<std::thread> threads;
    threads.reserve(kNumMethods);
    for (auto& runner : runners)
        threads.emplace_back([&runner]{ runner.run(); });

    // Ждём всех
    for (auto& t : threads)
        t.join();

    // Собираем результаты
    QVector<MethodResult> results;
    results.reserve(kNumMethods);
    for (const auto& r : runners)
        results.push_back(r.result());

    bestIdx = MethodSelector::selectBest(results);
    return results;
}

void ParallelRunner::runAllAsync(const SharedInputData& input)
{
    // Копируем данные — Approximator2D может потом изменится
    SharedInputData inputCopy = input;

    // Координирующий поток: запускает 4 рабочих, ждёт, испускает сигнал
    std::thread([this, inputCopy = std::move(inputCopy)]() mutable {
        int bestIdx = -1;
        QVector<MethodResult> results = this->runAll(inputCopy, bestIdx);

        // Вывод в консоль (из любого потока)
        MethodSelector::printResults(results, bestIdx);
        
        // emit из non-main thread: Qt сам обработает как надо и доставит
        emit resultsReady(results, bestIdx);
    }).detach();
}
