#ifndef AS_2D_PROCESS_CLASS_H
#define AS_2D_PROCESS_CLASS_H

#include<QVector>
#include<complex.h>

class AS_2D_process_class
{

private:
//    double sc;
//    double sf;
    int masksize;

//    QVector<QVector<double>> fp_kernel;
//    QVector<QVector<double>> q1_x_kernel;
//    QVector<QVector<double>> q1_y_kernel;
//    QVector<QVector<double>> q2_xx_kernel;
//    QVector<QVector<double>> q2_xy_kernel;
//    QVector<QVector<double>> q2_yy_kernel;
//    //
//    QVector<QVector<double>> f_s;
//    QVector<QVector<double>> f_plus;
//    QVector<QVector<double>> f_plus_minus;

public:
    AS_2D_process_class(int masksize);
    ~AS_2D_process_class();

    void fx_2DHilbertKernel(double s,
                            QVector<QVector<double>>& kernel_1st,
                            QVector<QVector<double>>& kernel_2nd);


    // Function to perform 2D FFT using FFTW
    QVector<QVector<std::complex<double>>> applyFFT2D(QVector<QVector<double>>& data);
    QVector<QVector<double>> applyIFFT2D(QVector<QVector<std::complex<double>>>& data);

    QVector<QVector<double>> convolve2DFreq(const QVector<QVector<double>>& image,
                                            const QVector<QVector<double>>& kernel);

    void AnalyticSignal_2dfft(QVector<QVector<double>>& Orientation,
                              QVector<QVector<double>>& phase,
                              QVector<QVector<double>>& amplitude,
                              QVector<QVector<double>>& apexangle,
                              const QVector<QVector<double>>& img);

    QVector<QVector<double>> y_mesh;
    QVector<QVector<double>> x_mesh;

};

#endif // AS_2D_PROCESS_CLASS_H
