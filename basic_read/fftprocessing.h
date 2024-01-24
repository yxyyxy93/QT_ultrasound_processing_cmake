#ifndef FFTPROCESSING_H
#define FFTPROCESSING_H

#include <QVector>
#include <complex>

class FFTProcessing {
public:
    static QVector<std::complex<double>> processFFT(QVector<double> tempVector);
};

#endif // FFTPROCESSING_H
