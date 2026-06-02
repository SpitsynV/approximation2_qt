#include "method_selector.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>

// ─────────────────────────────────────────────────────────────────────────────

int MethodSelector::selectBest(QVector<MethodResult>& results)
{
    if (results.isEmpty())
        return -1;

    double maxErr  = std::numeric_limits<double>::min();
    double maxTime = std::numeric_limits<double>::min();

    for (const auto& r : results) {
        if (!r.valid) continue;
        maxErr  = std::max(maxErr,  r.maxError);
        maxTime = std::max(maxTime, r.elapsedMs);
    }

    // Защита
    if (maxErr  < 1e-300) maxErr  = 1.0;
    if (maxTime < 1e-300) maxTime = 1.0;

    // Вычисляем score для каждого метода
    for (auto& r : results) {
        if (!r.valid) {
            r.score = std::numeric_limits<double>::max();
            continue;
        }
        const double normErr  = r.maxError  / maxErr;
        const double normTime = r.elapsedMs / maxTime;
        r.score = kAlpha * normErr + (1.0 - kAlpha) * normTime;
    }

    // Находим индекс лучшего метода (минимальный score)
    int bestIdx = -1;
    double bestScore = std::numeric_limits<double>::max();
    for (int i = 0; i < results.size(); ++i) {
        if (results[i].valid && results[i].score < bestScore) {
            bestScore = results[i].score;
            bestIdx   = i;
        }
    }
    return bestIdx;
}
// ─────────────────────────────────────────────────────────────────────────────
//Эту штуку можно запихать в QWidget
QString MethodSelector::toHtmlTable(const QVector<MethodResult>& results, int bestIdx)
{
    QString html;
    html += QStringLiteral(
        "<style>"
        "table{border-collapse:collapse;font-size:12px;font-family:monospace}"
        "th{background:#2d5a8e;color:white;padding:4px 8px;text-align:left}"
        "td{padding:3px 8px;border-bottom:1px solid #ccc}"
        ".best{background:#d4f0c0;font-weight:bold}"
        "</style>"
        "<table>"
        "<tr><th>Метод</th><th>maxError</th><th>Время (мс)</th><th>Score</th></tr>"
    );

    for (int i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        const bool  isBest = (i == bestIdx);
        const QString rowClass = isBest ? QStringLiteral(" class=\"best\"") : QString();

        html += QStringLiteral("<tr%1><td>%2%3</td><td>%4</td><td>%5</td><td>%6</td></tr>")
                    .arg(rowClass)
                    .arg(isBest ? QStringLiteral("(!) ") : QStringLiteral(""))
                    .arg(r.name.toHtmlEscaped())
                    .arg(r.maxError,  0, 'e', 3)
                    .arg(r.elapsedMs, 0, 'f', 2)
                    .arg(r.score,     0, 'f', 4);
    }

    html += QStringLiteral("</table>");

    if (bestIdx >= 0 && bestIdx < results.size()) {
        html += QStringLiteral(
            "<p style='margin:4px 0 0 0;font-size:12px;font-family:monospace;"
            "color:#1a6e1a;font-weight:bold'>(!) Лучший: %1</p>"
        ).arg(results[bestIdx].name.toHtmlEscaped());
    }

    return html;
}




// ─────────────────────────────────────────────────────────────────────────────

void MethodSelector::printResults(const QVector<MethodResult>& results, int bestIdx)
{
  using namespace std;

    const std::string sep(72, '-');

    cout << "\n" << sep << "\n";
    cout << "  ПАРАЛЛЕЛЬНОЕ СРАВНЕНИЕ МЕТОДОВ ИНТЕРПОЛЯЦИИ\n";
    cout << "  Критерий: score = " << kAlpha << "*E/Emax + "
         << (1.0 - kAlpha) << "*T/Tmax  (меньше = лучше)\n";
    cout << sep << "\n";
    cout << left  << setw(28) << "  Метод"
         << right << setw(16) << "maxError"
         << setw(12) << "T (мс)"
         << setw(10) << "Score"
         << "\n";
    cout << sep << "\n";

    for (int i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        std::string mark = (i == bestIdx) ? " (!)" : "  ";
        cout << left  << setw(28) << (mark + r.name.toStdString())
             << right << setw(16) << fixed << setprecision(6) << r.maxError
             << setw(12) << setprecision(3)  << r.elapsedMs
             << setw(10) << setprecision(4)  << r.score
             << "\n";
    }

    cout << sep << "\n";
    if (bestIdx >= 0 && bestIdx < results.size()) {
        const auto& best = results[bestIdx];
        cout << "  (!)ЛУЧШИЙ МЕТОД: " << best.name.toStdString()
             << "  (error=" << std::scientific << setprecision(6) << best.maxError
             << ",  time=" << setprecision(3) << best.elapsedMs << " мс"
             << ",  score=" << setprecision(4) << best.score << ")\n";
    }
    cout << sep << "\n\n";
}

// ─────────────────────────────────────────────────────────────────────────────
