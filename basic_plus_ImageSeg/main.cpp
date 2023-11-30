#include "ultrasound_cscan_seg.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);

    int fs = 250e6;
    double fx = 2e-3;
    double fy = 2e-3;
    QString default_fn = "default filename";
    ultrasound_cscan_seg process(nullptr,
                                 default_fn,
                                 fs,
                                 fx,
                                 fy);
    process.raise();
    process.activateWindow();
    process.show();
    return a.exec();
}
