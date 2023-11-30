#ifndef UTILS_H
#define UTILS_H

#endif // UTILS_H

#include <complex>
#include <QVector>
#include <qdebug.h>
#include <QLayout>
#include <random>

typedef QVector<double> Array1D;
typedef QVector<Array1D> Array2D;
typedef QVector<Array2D> Array3D;

QVector<std::complex<double>> analyticSignal(const QVector<double>& signal);

void fft(QVector<std::complex<double>>& signal);

void ifft(QVector<std::complex<double>>& signal);

// surface align utils
template <class T>
void shiftVector_1D(QVector<T>& signal, int k){
    for (int i=0; i<k; ++i){
        signal.append(signal.takeFirst());
    }
}

// 2d analytic-signal
// The first order 2D convolution kernels will be calculated by
double Kernel1(double x, double y, double s);

double Kernel2(double x, double y, double s);

void AnalyticSignal(double cx, double cy,
                    double& Orientation, double& Phase,
                    double& Amplitude, double& ApexAngle,
                    int n, double s_c, double s_f, double dx,
                    QVector<QVector<double>>& img);

QVector<double> min_max_2d(const QVector<QVector<double>> &data);

// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num);

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2);
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num);

// read files
QVector<double> readDatFile(const QString& filePath);

// trim dataset
template <class T>
QVector<QVector<QVector<T>>> trim3DData(const QVector<QVector<QVector<T>>>& data, int startI, int endI, int startJ, int endJ, int startK, int endK) {
    QVector<QVector<QVector<T>>> trimmedData;
    for (int i = startI; i <= endI; ++i) {
        QVector<QVector<T>> trimmed2DData;
        for (int j = startJ; j <= endJ; ++j) {
            QVector<T> trimmed1DData;
            for (int k = startK; k <= endK; ++k) {
                trimmed1DData.append(data[i][j][k]);
            }
            trimmed2DData.append(trimmed1DData);
        }
        trimmedData.append(trimmed2DData);
    }
    qDebug() << "data trimed";
    return trimmedData;
};
template<typename T>
QVector<T> downsampleVector(const QVector<T>& original, int rate)
{
    if (rate <= 0) {
        throw std::invalid_argument("Downsampling rate must be greater than 0");
    }

    QVector<T> downsampled;
    for (int i = 0; i < original.size(); i += rate) {
        downsampled.push_back(original[i]);
    }

    return downsampled;
}

// define operators for qvector of complex
QVector<QVector<QVector<double>>> abs(const QVector<QVector<QVector<std::complex<double>>>>& v);

// calculate statistics
QMap<int, int> calculateHistogram(const QVector<double>& data, int numBins);

// ******************** noises
// Function to generate random noise
double generateNoise(double amplitude);
// Function to add noise to a vector
void addNoiseToVector(std::vector<double>& vec, double amplitude);
double snrDbToLinear(double snrDb);
void addGaussianNoise(QVector<double>& signal, double snrDb);

template <typename T>
bool widgetExistsInLayout(QLayout* layout, const QString& propertyValue, const QString& propertyName = QString()) {
    if (layout) {
        for (int i = 0; i < layout->count(); ++i) {
            QWidget* widget = layout->itemAt(i)->widget();
            if (widget) {
                if (T* specificWidget = qobject_cast<T*>(widget)) {
                    // If propertyName is not provided, default to checking the objectName
                    QString valueToCheck = propertyName.isEmpty() ? specificWidget->objectName() : specificWidget->property(propertyName.toStdString().c_str()).toString();
                    if (valueToCheck == propertyValue) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
};

// ************** debug
void printWidgetInfo(QWidget *widget);
void printLayoutInfo(QLayout *layout);

// ********** fill nan
// Check if a value is NaN
bool isNaN(double value);
// Get the average of neighboring values for a given position in a 3D dataset
double getAverageOfNeighbors(const QVector<QVector<QVector<double>>>& data, int x, int y, int z);
// Fill NaN values in a 3D dataset with the average of their nearest neighbors
void fillNanWithNearestNeighbors(QVector<QVector<QVector<double>>>& data);
