#include "ultrasound_cscan_process.h"

#include <QApplication>
#include <QFileDialog>
#include <QFile>
#include <QTextStream>
#include <QDebug>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);

    // create the instance of the class
    int fs = 250e6;
    double fx = 2e-3;
    double fy = 2e-3;
    QString default_fn = "default filename";
    ultrasound_Cscan_process process(nullptr,
                                     default_fn,
                                     fs,
                                     fx,
                                     fy);
    process.show();
    return app.exec();
}
