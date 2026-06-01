#include <QApplication>
#include "plotwidget2d.h"
#include "approximator2d.h"
#include <iostream>
#include <cstdlib>
#include "method_result.h"
int main(int argc, char *argv[])
{
    if (argc < 10) {
        std::cerr << "Usage: " << argv[0]
                  << " a b c d nx ny mx my k" << std::endl;
        return EXIT_FAILURE;
    }

    double a  = std::stod(argv[1]);
    double b  = std::stod(argv[2]);
    double c  = std::stod(argv[3]);
    double d  = std::stod(argv[4]);
    int    nx = std::stoi(argv[5]);
    int    ny = std::stoi(argv[6]);
    int    mx = std::stoi(argv[7]);
    int    my = std::stoi(argv[8]);
    int    k  = std::stoi(argv[9]);

    QApplication app(argc, argv);
    qRegisterMetaType<MethodResult>("MethodResult");
    qRegisterMetaType<QVector<MethodResult>>("QVector<MethodResult>");

    Approximator2D approx(a, b, c, d, nx, ny, mx, my, k);

    PlotWidget2D plot(&approx);
    plot.setWindowTitle("Approx 2D ");
    plot.show();
    qRegisterMetaType<QVector<MethodResult>>();

    return app.exec();
}