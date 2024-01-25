#ifndef ULTRASOUND_CSCAN_SEG_H
#define ULTRASOUND_CSCAN_SEG_H

#include <QMainWindow>

#include <QtWidgets>

#include "../basic_read/ultrasound_cscan_process.h"
#include <memory>

#include "customgraphicsscene.h"

QT_BEGIN_NAMESPACE
namespace Ui { class ultrasound_cscan_seg; }
QT_END_NAMESPACE

class ultrasound_cscan_seg : public ultrasound_Cscan_process
{
    Q_OBJECT

public:
    ultrasound_cscan_seg(QWidget *parent,
                         QString fn,
                         int fs,
                         double fx,
                         double fy);
    ~ultrasound_cscan_seg();

    // functions for image manipulation
    template <typename T>
    QImage convertToImage(const QVector<QVector<QVector<T>>>& data, int z);

private slots:
    //regular
    void handleButton_Cscan();

    //  Drawing and Recording Arbitrary Lines:
    void handleButton_draw();
    void onLineDrawn(const QLineF &line);
    void saveImage();

    // slot on page 3
    void handleButton_multiSNR();
    void handleButtonSegmentDataset();
    void segmentAndSaveData(const QVector<QVector<QVector<std::complex<double>>>> data3d);
    void selectFolder();

    // slot on page 4
    void handleCalculateFFT();
    void handleCalculateCepstrum();

public slots:
    void closeDrawnArea();

private:
//    Ui::ultrasound_cscan_seg *ui;
    QVBoxLayout *layout2;
    QVBoxLayout *layout3;

    // Qt plots
    QCustomPlot *customPlot1_page2;
    QCustomPlot *customPlot2_page2;
    QScrollBar *scrollBarZ_page2;

    // images properties
    QVector<QVector<double>> img_Cscan;
    QPointF startPoint;
    CustomGraphicsScene *scene=nullptr;
    QGraphicsView *embeddedView;
    QVector<QVector<QVector<bool>>> C_scan_mask;

    // basic propertes
    int x_size;
    int y_size;
    int z_size;

    //  functions on page 3
    QProgressBar* progressBarPage3;
    QLineEdit* multisnrInput;
    QLineEdit* downsampleRateInput;
    QLineEdit* cropSignalInput;
    QComboBox* myComboBox_savepattern;

    QLabel *folderLabel;
    void processFolder(const QString &path);
    void readCropValues(int& start, int& end);

    // functions on page 4
    QVector<QVector<QVector<std::complex<double>>>> fft3d;
    QVector<QVector<QVector<std::complex<double>>>> cepstrum3d;
};
#endif // ULTRASOUND_CSCAN_SEG_H





