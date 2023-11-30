#include "utils.h"

#include <cmath>
#include <complex>
#include <QVector>
#include "fftw3.h"

#include <QDir>
#include <QFile>
#include <QDataStream>
#include <QDebug>
#include <QPushButton>
#include <QLayout>
#include <QDebug>

#include <algorithm> // for std::count_if
#include <functional> // for std::bind, std::placeholders, and std::less

#include <limits>
#include <cmath>  // for std::isnan

// analytic-signal calculation
QVector<std::complex<double>> analyticSignal(const QVector<double>& signal)
{
    int N = signal.size();

    // Compute the FFT of the signal
    QVector<std::complex<double>> signalFFT(N);
    for (int i = 0; i < N; i++) {
        signalFFT[i] = std::complex<double>(signal[i], 0.0);
    }
    fft(signalFFT); // Replace with your own FFT function

    // Compute the frequency domain representation of the Hilbert transform
    QVector<std::complex<double>> hilbertFFT(N);
    hilbertFFT[0] = std::complex<double>(0.0, 0.0);
    for (int i = 1; i < N/2; i++) {
        hilbertFFT[i] = std::complex<double>(0.0, -1.0/(M_PI*i));
        hilbertFFT[N-i] = std::conj(hilbertFFT[i]);
    }
    hilbertFFT[N/2] = std::complex<double>(0.0, 0.0);

    // Multiply the FFT of the signal by the frequency domain representation of the Hilbert transform
    QVector<std::complex<double>> analyticFFT(N);
    for (int i = 0; i < N/2; i++) {
        analyticFFT[i] = signalFFT[i] * std::complex<double>(2.0, 0.0);
    }
    for (int i = N/2; i < N; i++) {
        analyticFFT[i] = std::complex<double>(0.0, 0.0);
    }

    // Compute the inverse FFT of the analytic signal
    ifft(analyticFFT); // Replace with your own IFFT function

    return analyticFFT;
}

void fft(QVector<std::complex<double>>& signal)
{
    int N = signal.size();

    // Create a FFTW plan for forward FFT
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&signal[0]),
            reinterpret_cast<fftw_complex*>(&signal[0]),
            FFTW_FORWARD, FFTW_ESTIMATE);

    // Execute the FFTW plan
    fftw_execute(plan);

    // Destroy the FFTW plan
    fftw_destroy_plan(plan);
}

void ifft(QVector<std::complex<double>>& signal)
{
    int N = signal.size();
    // Create a FFTW plan for backward FFT
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&signal[0]),
            reinterpret_cast<fftw_complex*>(&signal[0]),
            FFTW_BACKWARD, FFTW_ESTIMATE);

    // Execute the FFTW plan
    fftw_execute(plan);

    // Normalize the output of IFFT
    for (int i = 0; i < N; i++) {
        signal[i] /= N;
    }

    // Destroy the FFTW plan
    fftw_destroy_plan(plan);
}

// 2d analytic-signal
// The first order 2D convolution kernels will be calculated by
double Kernel1(double x, double y, double s)
{
    double ss = pow(s, 2);
    double kk = pow(x, 2) + pow(y, 2);
    return 1 / (2 * M_PI * pow(ss + kk, 1.5));
}
//the second order 2D convolution kernels in spatial domain(x, y) with scale space parameter s will be determined by

double Kernel2(double x, double y, double s)
{
    double ss = pow(s, 2);
    double kk = pow(x, 2) + pow(y, 2);
    double d = pow(kk, 2) * pow(ss + kk, 1.5) * 2 * M_PI;
    return (d == 0) ? 0 : (s * (2 * ss + 3 * kk) - 2 * pow(ss + kk, 1.5)) / d;
}

void AnalyticSignal(double cx, double cy,
                    double& Orientation, double& Phase,
                    double& Amplitude, double& ApexAngle,
                    int n, double s_c, double s_f, double dx,
                    QVector<QVector<double>>& img)
{
    double f_p = 0, f_x = 0, f_y = 0, f_xx = 0, f_xy = 0, f_yy = 0;
    //2D convolution
    for (int x = -n;x <= n;x += dx)
        for (int y = -n;y <= n;y += dx)
        {
            double t = img[x + cx][y + cy] * pow(dx, 2);
            double pf = t * Kernel1(x, y, s_f);
            double pc = t * Kernel1(x, y, s_c);
            double k = t * (Kernel2(x, y, s_f) - Kernel2(x, y, s_c));
            //signal in Poisson scale space
            f_p += s_f * pf - s_c * pc;
            //first order Hilbert transform
            f_x += x * (pf - pc);
            f_y += y * (pf - pc);
            //second order Hilbert transform
            f_xx += x * x * k;
            f_yy += y * y * k;
            f_xy += x * y * k;
        }
    double f_pm = 0.5 * (f_xx - f_yy);
    double f_s = 0.5 * f_p;
    double f_plus = f_xy;
    double e = sqrt(pow(f_pm, 2) + pow(f_plus, 2)) / fabs(f_s);
    double q = (pow(f_x, 2) + pow(f_y, 2)) * 2 / (1 + e);
    Phase = atan2(sqrt(q), f_p);
    Orientation = ((int)Phase == 0)
            ? 0.5 * atan2(f_plus, f_pm) + M_PI / 2
            : atan2(f_y, f_x);
    Amplitude = 0.5 * sqrt(pow(f_p, 2) + q);
    ApexAngle = atan2(sqrt(pow(f_s, 2) - pow(f_plus, 2) -
                           pow(f_pm, 2)), sqrt(pow(f_plus, 2) + pow(f_pm, 2)));
}

QVector<double> min_max_2d(const QVector<QVector<double>> &data){
    // initialize the minimum and maximum values to the first element in the vector
    double minVal = data[0][0];
    double maxVal = data[0][0];

    // iterate over all the elements in the vector and update the minimum and maximum values
    for (int i = 0; i < data.size(); ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            if (data[i][j] < minVal) {
                minVal = data[i][j];
            }
            if (data[i][j] > maxVal) {
                maxVal = data[i][j];
            }
        }
    }
    QVector<double> results;
    results.push_back(minVal);
    results.push_back(maxVal);

    return results;
}

// To overload the +, -, * operators with the input type QVector<QVector<double>>, point-wise operation
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator+(const QVector<QVector<double>>& v1, const double& num){
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] + num;
        }
    }

    return result;
}

QVector<QVector<double>> operator-(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] - v2[i][j];
        }
    }

    return result;
}

QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const QVector<QVector<double>>& v2) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * v2[i][j];
        }
    }

    return result;
}
QVector<QVector<double>> operator*(const QVector<QVector<double>>& v1, const double& num) {
    int numRows = v1.size();
    int numCols = v1[0].size();

    QVector<QVector<double>> result(numRows, QVector<double>(numCols));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            result[i][j] = v1[i][j] * num;
        }
    }

    return result;
}

// read files
QVector<double> readDatFile(const QString& filePath) {
    QFile file(filePath);
    if (!file.open(QIODevice::ReadOnly)) {
        qWarning() << "Cannot open file for reading:" << file.errorString();
        return QVector<double>();
    }

    QDataStream in(&file);
    // Set the byte order to LittleEndian if the file was written on a little-endian machine
    in.setByteOrder(QDataStream::LittleEndian);

    QVector<double> data;
    while (!in.atEnd()) {
        double Value;
        in >> Value;
        data.append(Value);
    }

    return data;
}

// define functions for vector of complex
QVector<QVector<QVector<double>>> abs(const QVector<QVector<QVector<std::complex<double>>>>& v){
    int numRows = v.size();
    int numCols = v[0].size();
    int numdeps = v[0][0].size();

    QVector<QVector<QVector<double>>> result(numRows, QVector<QVector<double>>(numCols, QVector<double>(numdeps)));

    for (int i = 0; i < numRows; i++) {
        for (int j = 0; j < numCols; j++) {
            for (int q = 0; q < numdeps; q++){
                result[i][j][q] = sqrt(pow(v[i][j][q].real(), 2) + pow(v[i][j][q].imag(), 2));
            }
        }
    }

    return result;
}

// calculate statistics
QMap<int, int> calculateHistogram(const QVector<double>& data, int numBins)
{
    QMap<int, int> histogram;

    double min = *std::min_element(data.begin(), data.end());
    double max = *std::max_element(data.begin(), data.end());
    double binSize = (max - min) / numBins;

    for (double value : data)
    {
        int bin = (value - min) / binSize;
        histogram[bin]++;
    }

    return histogram;
}

// ******************** noises
// Function to generate random noise
double generateNoise(double amplitude) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-amplitude, amplitude);
    return dis(gen);
}
// Function to add noise to a vector
void addNoiseToVector(std::vector<double>& vec, double amplitude) {
    for (auto& element : vec) {
        element += generateNoise(amplitude);
    }
}
double snrDbToLinear(double snrDb) {
    return std::pow(10.0, snrDb / 10.0);
}
void addGaussianNoise(QVector<double>& signal, double snrDb) {
    double snrLinear = snrDbToLinear(snrDb);
    double signalPower = std::accumulate(signal.begin(), signal.end(), 0.0,
                                         [](double sum, double value) { return sum + value * value; }) / signal.size();
    double noisePower = signalPower / snrLinear;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0, std::sqrt(noisePower));

    for (auto& value : signal) {
        value += d(gen);
    }
}

// ************** debug
void printWidgetInfo(QWidget *widget) {
    if (widget) {
        qDebug() << "Widget:" << widget;
        qDebug() << "Size Policy:" << widget->sizePolicy().horizontalPolicy()
                 << widget->sizePolicy().verticalPolicy();
        qDebug() << "Minimum Size:" << widget->minimumSize();
        qDebug() << "Maximum Size:" << widget->maximumSize();
        qDebug() << "Actual Size:" << widget->size();
        qDebug() << "--------------------------------";
    }
}
void printLayoutInfo(QLayout *layout) {
    if (layout) {
        qDebug() << "Layout:" << layout;
        qDebug() << "Number of items:" << layout->count();

        for (int i = 0; i < layout->count(); ++i) {
            QLayoutItem *item = layout->itemAt(i);
            if (item->widget()) {
                qDebug() << "Item" << i << "is a widget.";
                printWidgetInfo(item->widget());
            } else if (item->spacerItem()) {
                qDebug() << "Item" << i << "is a spacer item.";
                qDebug() << "Size:" << item->spacerItem()->sizeHint();
            }
            // Add more checks here if you have other types of items
        }
        qDebug() << "================================";
    }
}

// **************** fill in nan
bool isNaN(double value) {
    return std::isnan(value);
}
double getAverageOfNeighbors(const QVector<QVector<QVector<double>>>& data, int x, int y, int z) {
    double sum = 0.0;
    int count = 0;
    int dx[] = {-1, 0, 1, 0, 0, 0};
    int dy[] = {0, -1, 0, 1, 0, 0};
    int dz[] = {0, 0, 0, 0, -1, 1};

    for (int i = 0; i < 6; ++i) {
        int nx = x + dx[i];
        int ny = y + dy[i];
        int nz = z + dz[i];

        if (nx >= 0 && nx < data.size() &&
            ny >= 0 && ny < data[nx].size() &&
            nz >= 0 && nz < data[nx][ny].size() &&
            !isNaN(data[nx][ny][nz])) {
            sum += data[nx][ny][nz];
            ++count;
        }
    }

    if (count > 0) {
        return sum / count;
    } else {
        return 0;  // Default value if no neighbors are available
    }
}
void fillNanWithNearestNeighbors(QVector<QVector<QVector<double>>>& data) {
    for (int x = 0; x < data.size(); ++x) {
        for (int y = 0; y < data[x].size(); ++y) {
            for (int z = 0; z < data[x][y].size(); ++z) {
                if (isNaN(data[x][y][z])) {
                    qDebug() << "find Nan";
                    data[x][y][z] = getAverageOfNeighbors(data, x, y, z);
                }
            }
        }
    }
}

