#include "ultrasound_Cscan_process_2das.h"

#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    // create the instance of the class
    int fs = 200e6;
    double fx = 2.5e-3;
    double fy = 2.5e-3;
    QString default_fn = "default filename";
    ultrasound_cscan_process_2DAS process(nullptr,
                                          default_fn,
                                          fs,
                                          fx,
                                          fy);
    process.show();
    return a.exec();
}
