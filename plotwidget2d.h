#pragma once
#include <QWidget>
#include <QVector>
#include <QLabel>
#include <QPushButton>
#include "approximator2d.h"
#include "parallel_runner.h"
#include "method_result.h"
#include <QLabel>
#include <QPushButton>

class PlotWidget2D : public QWidget
{
    Q_OBJECT
public:
    explicit PlotWidget2D(Approximator2D *approx, QWidget *parent = nullptr);
    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;

protected:
    void paintEvent(QPaintEvent *) override;
    void keyPressEvent(QKeyEvent *event) override;

    void resizeEvent(QResizeEvent *e);

private:
    Approximator2D *m_approx;
    double m_lastMaxAbs;

    // Проецируем точку (xw, yw, zw) в экранные координаты
    // с учётом угла поворота m_approx->angle() вокруг OZ
    QPointF project(double xw, double yw, double zw,
                    double cx, double cy,
                    double scaleXY, double scaleZ,
                    double xmid, double ymid, double zmid) const;
    


    //

    // ── Параллельный запуск ──────────────────────────────────────
    ParallelRunner* m_runner       = nullptr;
    QLabel*         m_resultsLabel = nullptr;  // панель результатов
    QPushButton*    m_runBtn       = nullptr;  // кнопка "Сравнить методы"
    private slots:
    void onRunParallel();
    void onResultsReady(QVector<MethodResult> results, int bestIdx);
};