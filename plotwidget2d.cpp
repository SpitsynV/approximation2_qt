#include "plotwidget2d.h"
#include <QPainter>
#include <QKeyEvent>
#include <cmath>
#include <cstdio>
#include <algorithm>

PlotWidget2D::PlotWidget2D(Approximator2D *approx, QWidget *parent)
    : QWidget(parent), m_approx(approx), m_lastMaxAbs(-1.0)
{
    setFocusPolicy(Qt::StrongFocus);
}

QSize PlotWidget2D::minimumSizeHint() const { return QSize(400, 300); }
QSize PlotWidget2D::sizeHint()        const { return QSize(800, 600); }

// Проекция мирового (xw,yw,zw) → экранное (u,v).
// Вращение вокруг OZ на угол theta, затем косоугольная проекция:
//   xr = (xw-xmid)*cos(θ) - (yw-ymid)*sin(θ)
//   yr = (xw-xmid)*sin(θ) + (yw-ymid)*cos(θ)
//   u  = cx + scaleXY * xr - 0.5*scaleXY * yr
//   v  = cy - scaleZ  * (zw-zmid) - 0.4*scaleXY * yr
QPointF PlotWidget2D::project(double xw, double yw, double zw,
                              double cx, double cy,
                              double scaleXY, double scaleZ,
                              double xmid, double ymid, double zmid) const
{
    double theta = m_approx->angle();
    double dx = xw - xmid;
    double dy = yw - ymid;
    double xr =  dx * std::cos(theta) - dy * std::sin(theta);
    double yr =  dx * std::sin(theta) + dy * std::cos(theta);
    double u = cx + scaleXY * xr - 0.5  * scaleXY * yr;
    double v = cy - scaleZ  * (zw - zmid) - 0.4 * scaleXY * yr;
    return {u, v};
}

void PlotWidget2D::paintEvent(QPaintEvent *)
{
    if (!m_approx) return;
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing, false);  // быстрее без AA

    auto plotFunc = m_approx->getPlotFunc();

    // ── Текущая область (с учётом масштаба s) ─────────────────────────
    const int s = m_approx->scaleExp();
    const double factor = 1.0 / (1 << s);

    double xmid = (m_approx->a() + m_approx->b()) * 0.5;
    double ymid = (m_approx->c() + m_approx->d()) * 0.5;
    double xhalf = (m_approx->b() - m_approx->a()) * 0.5 * factor;
    double yhalf = (m_approx->d() - m_approx->c()) * 0.5 * factor;

    double xMin = xmid - xhalf, xMax = xmid + xhalf;
    double yMin = ymid - yhalf, yMax = ymid + yhalf;

    // ── Визуализационная сетка mx×my ──────────────────────────────────
    // ограничиваем отображаемое число ячеек размером экрана
    const int mx = m_approx->mx();
    const int my = m_approx->my();
    const int dispX = std::min(mx, std::max(2, width()  / 4));
    const int dispY = std::min(my, std::max(2, height() / 4));

    std::vector<double> gx(dispX), gy(dispY);
    for (int i = 0; i < dispX; i++)
        gx[i] = xMin + i * (xMax - xMin) / (dispX - 1);
    for (int j = 0; j < dispY; j++)
        gy[j] = yMin + j * (yMax - yMin) / (dispY - 1);

    // вычисляем значения на сетке
    std::vector<double> fgrid(dispX * dispY);
    double Fmin = 0.0, Fmax = 0.0;
    bool firstVal = true;
    for (int i = 0; i < dispX; i++) {
        for (int j = 0; j < dispY; j++) {
            double v = plotFunc(gx[i], gy[j]);
            fgrid[i * dispY + j] = v;
            if (std::isfinite(v)) {
                if (firstVal) { Fmin = Fmax = v; firstVal = false; }
                else { Fmin = std::min(Fmin, v); Fmax = std::max(Fmax, v); }
            }
        }
    }

    double maxAbs = std::max(std::fabs(Fmin), std::fabs(Fmax));
    if (std::fabs(maxAbs - m_lastMaxAbs) > 1e-16 || m_lastMaxAbs < 0) {
        fprintf(stderr, "max{|Fmin|,|Fmax|} = %g  (s=%d)\n", maxAbs, s);
        m_lastMaxAbs = maxAbs;
    }

    double zmid = 0.5 * (Fmin + Fmax);
    double zrange = Fmax - Fmin;
    if (zrange < 1e-14) zrange = 1.0;

    // ── Масштабирование проекции ───────────────────────────────────────
    // подбираем scaleXY и scaleZ так, чтобы поверхность вписалась в окно
    double xyrange = std::max(xMax - xMin, yMax - yMin);
    double scaleXY = std::min(width(), height()) * 0.4 / (xyrange > 0 ? xyrange : 1.0);
    double scaleZ  = std::min(width(), height()) * 0.35 / zrange;

    double cx = width()  * 0.5;
    double cy = height() * 0.55;

    // ── Каркасная отрисовка поверхности ───────────────────────────────
    painter.setPen(QPen(Qt::blue, 1));

    // линии вдоль x (при фиксированном j)
    for (int j = 0; j < dispY; j++) {
        QPointF prev = project(gx[0], gy[j], fgrid[0 * dispY + j],
                               cx, cy, scaleXY, scaleZ, xmid, ymid, zmid);
        for (int i = 1; i < dispX; i++) {
            double v = fgrid[i * dispY + j];
            if (!std::isfinite(v)) { prev = project(gx[i], gy[j], 0, cx, cy, scaleXY, scaleZ, xmid, ymid, zmid); continue; }
            QPointF cur = project(gx[i], gy[j], v, cx, cy, scaleXY, scaleZ, xmid, ymid, zmid);
            painter.drawLine(prev, cur);
            prev = cur;
        }
    }

    // линии вдоль y (при фиксированном i)
    painter.setPen(QPen(QColor(0, 150, 200), 1));
    for (int i = 0; i < dispX; i++) {
        QPointF prev = project(gx[i], gy[0], fgrid[i * dispY + 0],
                               cx, cy, scaleXY, scaleZ, xmid, ymid, zmid);
        for (int j = 1; j < dispY; j++) {
            double v = fgrid[i * dispY + j];
            if (!std::isfinite(v)) { prev = project(gx[i], gy[j], 0, cx, cy, scaleXY, scaleZ, xmid, ymid, zmid); continue; }
            QPointF cur = project(gx[i], gy[j], v, cx, cy, scaleXY, scaleZ, xmid, ymid, zmid);
            painter.drawLine(prev, cur);
            prev = cur;
        }
    }

    // ── Информационная строка ──────────────────────────────────────────
    painter.setPen(Qt::black);
    QFont font = painter.font();
    font.setPointSize(10);
    painter.setFont(font);

    double angleDeg = m_approx->angle() * 180.0 / M_PI;
    QString info = QString("k=%1 %2  nx=%3 ny=%4  s=%5  p=%6  angle=%7°  %8  max|F|=%9")
        .arg(m_approx->k())
        .arg(m_approx->functionName())
        .arg(m_approx->nx())
        .arg(m_approx->ny())
        .arg(s)
        .arg(m_approx->p())
        .arg(angleDeg, 0, 'f', 0)
        .arg(m_approx->getPlotName())
        .arg(maxAbs, 0, 'g', 4);
    painter.drawText(10, 20, info);
}

void PlotWidget2D::keyPressEvent(QKeyEvent *event)
{
    if (!m_approx) { QWidget::keyPressEvent(event); return; }

    bool needRebuild = false;

    switch (event->key()) {
        //debug
        case Qt::Key_D:
    {
    //const auto &x = m_approx->getX();
    //const auto &f = m_approx->getF();
    //fprintf(stderr, "=== Debug: n=%d, a=%.6f, b=%.6f ===\n", m_approx->n(), m_approx->a(), m_approx->b());
    
    //fprintf(stderr, "Max error method 1 = %e\n", m_approx->getMaxError1());
    //fprintf(stderr, "Max error method 2 = %e\n", m_approx->getMaxError2());
    // не меняем состояние и не перерисовываем
    return;
    }



    // 0 — циклически менять k (номер функции)
    case Qt::Key_0:
        m_approx->nextK();
        needRebuild = true;
        break;

    // 1 — циклически менять состав графиков
    case Qt::Key_1:
        m_approx->nextGraphMode();
        break;

    // 2 — увеличить масштаб (растяжение XY в 2 раза)
    case Qt::Key_2:
        m_approx->setScaleExp(m_approx->scaleExp() + 1);
        break;

    // 3 — уменьшить масштаб
    case Qt::Key_3:
        if (m_approx->scaleExp() > 0)
            m_approx->setScaleExp(m_approx->scaleExp() - 1);
        break;

    // 4 — увеличить nx,ny в 2 раза
    case Qt::Key_4:
        m_approx->setNx(m_approx->nx() * 2);
        m_approx->setNy(m_approx->ny() * 2);
        needRebuild = true;
        break;

    // 5 — уменьшить nx,ny в 2 раза
    case Qt::Key_5:
        m_approx->setNx(m_approx->nx() / 2);
        m_approx->setNy(m_approx->ny() / 2);
        needRebuild = true;
        break;

    // 6 — увеличить возмущение p на 1
    case Qt::Key_6:
        m_approx->setP(m_approx->p() + 1);
        needRebuild = true;
        break;

    // 7 — уменьшить возмущение p на 1
    case Qt::Key_7:
        if (m_approx->p() > 0) m_approx->setP(m_approx->p() - 1);
        else                    m_approx->setP(m_approx->p() - 1);
        needRebuild = true;
        break;

    // 8 — вращать по часовой стрелке вокруг OZ на π/12
    case Qt::Key_8:
        m_approx->setAngle(m_approx->angle() + M_PI / 12.0);
        break;

    // 9 — вращать против часовой стрелки
    case Qt::Key_9:
        m_approx->setAngle(m_approx->angle() - M_PI / 12.0);
        break;

    default:
        QWidget::keyPressEvent(event);
        return;
    }

    if (needRebuild) {
        m_approx->rebuild();
        fprintf(stderr, "nx=%d ny=%d  p=%d\n",
                m_approx->nx(), m_approx->ny(), m_approx->p());
    }

    m_lastMaxAbs = -1.0;  // принудительно выводим max|F| при следующей отрисовке
    update();
}