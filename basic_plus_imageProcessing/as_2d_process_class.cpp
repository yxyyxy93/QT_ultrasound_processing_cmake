#include "as_2d_process_class.h"
#include <cmath>
#include <QDebug>
#include <fftw3.h>
#include <complex.h>

AS_2D_process_class::AS_2D_process_class(int masksize):
    masksize(masksize) {
    // define the meshgrids
    // the real size is 2*masksize+1
    QVector<double> oneline_x;
    QVector<double> oneline_y;
    for (int i = -masksize; i <= masksize; ++i) {
        oneline_x.clear();
        oneline_y.clear();
        for (int j = -masksize; j <= masksize; ++j) {
            oneline_x.push_back(static_cast<double>(i));
            oneline_y.push_back(static_cast<double>(j));
        }
        this->x_mesh.push_back(oneline_x);
        this->y_mesh.push_back(oneline_y);
    }
}

AS_2D_process_class::~AS_2D_process_class(){

}

void AS_2D_process_class::fx_2DHilbertKernel(double s,
                                             QVector<QVector<double>>& kernel_1st,
                                             QVector<QVector<double>>& kernel_2nd) {
    double ss = pow(s, 2);
    double kk;
    double d;
    QVector<double> kernel_1st_oneline;
    QVector<double> kernel_2nd_oneline;
    for (int i = -this->masksize; i <= this->masksize; ++i) {
        kernel_1st_oneline.clear();
        kernel_2nd_oneline.clear();
        for (int j = -this->masksize; j <= this->masksize; ++j) {
            kk = pow(i, 2) + pow(j, 2);
            // 1st kernel
            kernel_1st_oneline.push_back(1 / (2*M_PI*pow(ss+kk, 1.5)));
            // 2nd kernel
            d = pow(kk, 2) * pow(ss+kk, 1.5) * 2 * M_PI;
            kernel_2nd_oneline.push_back(d==0? 0:-(s * (2*ss + 3*kk) - 2*pow(ss+kk, 1.5)) / d);
        }
        kernel_1st.push_back(kernel_1st_oneline);
        kernel_2nd.push_back(kernel_2nd_oneline);
    }
}

QVector<QVector<std::complex<double>>> AS_2D_process_class::
applyFFT2D(QVector<QVector<double>>& data)
{
    QVector<QVector<std::complex<double>>> output(data.size(), QVector<std::complex<double>>(data[0].size()));
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_plan plan = fftw_plan_dft_2d(data.size(), data[0].size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    // copy input data to fftw_complex array
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            in[i*data[0].size()+j][0] = data[i][j];
            in[i*data[0].size()+j][1] = 0.0;
        }
    }
    fftw_execute(plan);
    // copy output data to QVector<QVector<std::complex<double>>> format
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            output[i][j] = std::complex<double>(out[i*data[0].size()+j][0], out[i*data[0].size()+j][1]);
        }
    }
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return output;
}

QVector<QVector<double>> AS_2D_process_class::
applyIFFT2D(QVector<QVector<std::complex<double>>>& data)
{
    QVector<QVector<double>> output(data.size(), QVector<double>(data[0].size()));
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * data.size() * data[0].size());
    fftw_plan plan = fftw_plan_dft_2d(data.size(), data[0].size(), in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    // copy input data to fftw_complex array
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            in[i*data[0].size()+j][0] = data[i][j].real();
            in[i*data[0].size()+j][1] = data[i][j].imag();
        }
    }
    fftw_execute(plan);
    // copy output data to QVector<QVector<double>> format
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0].size(); j++) {
            output[i][j] = out[i*data[0].size()+j][0] / (data.size()*data[0].size());
        }
    }
    fftw_destroy_plan(plan);
    fftw_free(in);
    fftw_free(out);
    return output;
}

// Compute the 2D convolution of an image and a kernel in the frequency domain
QVector<QVector<double>> AS_2D_process_class::
convolve2DFreq(const QVector<QVector<double>>& image, const QVector<QVector<double>>& kernel)
{
    int numRows = image.size();
    int numCols = image[0].size();
    int kernelNumRows = kernel.size();
    int kernelNumCols = kernel[0].size();
    // Pad the input image and the kernel to avoid circular convolution
    int paddedNumRows = numRows + kernelNumRows - 1;
    int paddedNumCols = numCols + kernelNumCols - 1;
    QVector<QVector<double>> paddedImage(paddedNumRows, QVector<double>(paddedNumCols, 0.0));
    QVector<QVector<double>> paddedKernel(paddedNumRows, QVector<double>(paddedNumCols, 0.0));
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            paddedImage[i][j] = image[i][j];
        }
    }
    for (int i = 0; i < kernelNumRows; i++) {
        for (int j = 0; j < kernelNumCols; j++) {
            paddedKernel[i][j] = kernel[i][j];
        }
    }
    // Compute the IDFT of the product
    QVector<QVector<std::complex<double>>> Image_fft = this->applyFFT2D(paddedImage);
    QVector<QVector<std::complex<double>>> Kernel_fft = this->applyFFT2D(paddedKernel);
    // Multiply the frequency-domain representations of the image and the kernel
    QVector<QVector<std::complex<double>>> fftResult(paddedNumRows,
                                                     QVector<std::complex<double>>(paddedNumCols, 0.0));
    for (int i = 0; i < paddedNumRows; i++) {
        for (int j = 0; j < paddedNumCols; j++) {
            fftResult[i][j] = Image_fft[i][j] * Kernel_fft[i][j];
        }
    }
    QVector<QVector<double>> Result = this->applyIFFT2D(fftResult);
    // Crop the result to remove the padding
    QVector<QVector<double>> croppedResult(numRows, QVector<double>(numCols, 0.0));
    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            croppedResult[i][j] = Result[i + kernelNumRows / 2][j + kernelNumCols / 2];
        }
    }
    return croppedResult;
}
