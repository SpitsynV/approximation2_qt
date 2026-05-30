#pragma once
#include <QWidget>
#include "approximator2d.h"

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

private:
    Approximator2D *m_approx;
    double m_lastMaxAbs;

    // Проецируем точку (xw, yw, zw) в экранные координаты
    // с учётом угла поворота m_approx->angle() вокруг OZ
    QPointF project(double xw, double yw, double zw,
                    double cx, double cy,
                    double scaleXY, double scaleZ,
                    double xmid, double ymid, double zmid) const;
};