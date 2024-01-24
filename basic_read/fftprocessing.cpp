#include "fftprocessing.h"
#include "utils.h"

QVector<std::complex<double>> FFTProcessing::processFFT(QVector<double> tempVector) {
    QVector<std::complex<double>> fftVector = applyFFT1D(tempVector);
    QVector<std::complex<double>> half_fftVector(fftVector.begin(), fftVector.begin() + fftVector.size() / 2);

    QVector<double> abs_half_fftVector;
    for (const auto& value : half_fftVector) {
        abs_half_fftVector.push_back(std::abs(value));
    }

    QVector<std::complex<double>> fft2Vector = applyFFT1D(abs_half_fftVector);
    QVector<std::complex<double>> half_fft2Vector(fft2Vector.begin(), fft2Vector.begin() + fft2Vector.size() / 2);

    return half_fft2Vector;
}
