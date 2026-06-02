#pragma once
#include "shared_input.h"
#include "method_result.h"
#include <functional>
/// Использование:
///   MethodRunner runner(2, inputData);
///   std::thread t([&runner]{ runner.run(); });
///   t.join();
///   auto res = runner.result();
class MethodRunner {
public:
    explicit MethodRunner(int methodId, const SharedInputData& input);
    void run();
    const MethodResult& result() const { return m_result; }

private:
    int                    m_methodId;
    const SharedInputData& m_input;    // только чтение
    MethodResult           m_result;

    /// Разрешение тестовой сетки для оценки погрешности.
    static constexpr int kMinTestN = 50;

    double computeMaxError(std::function<double(double, double)> approxFn, std::function<double(double, double)> exactFn) const;

    MethodResult runMethod1(); // Метод Бесселя            (task3)
    MethodResult runMethod2(); // Разделённые разности      (task31)
    MethodResult runMethod3(); // Эрмит                     (task5)
    MethodResult runMethod4(); // МНК по Чебышеву          (task6)
};
