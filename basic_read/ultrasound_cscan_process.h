#ifndef ULTRASOUND_CSCAN_PROCESS_H
#define ULTRASOUND_CSCAN_PROCESS_H

#include <QPushButton>
#include <QtWidgets>
#include <QProgressBar>
#include <complex>
#include <cmath>

#include "qcustomplot.h"

class ultrasound_Cscan_process : public QWidget
{
    Q_OBJECT

public:
    // constructor
    ultrasound_Cscan_process(QWidget *parent = nullptr,
                             QString fn = "",
                             int fs = 250e6,
                             double fx = 2e-3,
                             double fy = 2e-3);
    ~ultrasound_Cscan_process();

    // visualization
    void setData();
    void updatePlot();

    friend class ultrasound_cscan_process_2DAS;
    friend class ultrasound_cscan_seg;

private slots:
    void handleButton_load();
    void handleButton_save();
    // void handleButton_loadraw();
    void handleButton_orthoslice();
    void handleButton_surface();
    void onCustomPlotClicked_Cscan(QMouseEvent* event);

    // in "handleButton_surface"
    void handleButton_plotsurface();
    void handleButton_alignsurface();

    // manipulate the dataset
    void handleButton_trim();
    void handleButton_addNoise();

    // delete dynamic allocations
    void clearAllDynamicMemory();

    // Cscan setting
    void updateCscanPlotSelection(int index);

    // manipulate widgets
    void addNewWidgetAndReorderLayout(QWidget* newWidget);

private:
    // Qt interfaces for load file
    QPushButton *myButton_load;
    QPushButton *myButton_save;
    void processFile(const QFileInfo &fileInfo);

    // QPushButton *myButton_loadraw;
    QPushButton *myButton_orthoslice;
    QPushButton *myButton_surface;

    QPushButton* myButton_plotsurface = nullptr;
    QPushButton* myButton_alignsurface = nullptr;

    QPushButton* myButton_addNoise=nullptr;
    QLineEdit* snrInput=nullptr;

    QProgressBar *m_progressBar;
    // basic properties
    QString fn;
    int fs;
    int fx;
    int fy;

    QVector<QVector<double>> Front_surface_idx;
    QVector<QVector<double>> Front_surface_val;

    // Qt plots
    QCustomPlot *customPlot1;
    QCustomPlot *customPlot2;
    QCustomPlot *customPlot3;
    QScrollBar *scrollBarX;
    QScrollBar *scrollBarY;
    QScrollBar *scrollBarZ;

    QCustomPlot *customPlot_Ascan;

    QCustomPlot *customPlot_frontI;
    QCustomPlot *customPlot_frontV;

    // dynamic memory management
    QVector<QCustomPlot*> customPlots;
    QVector<QScrollBar*> scrollBars;
    QVector<QLabel*> labels;
    QVector<QWidget*> widgets;
    QVector<QHBoxLayout*> hLayouts;
    QVector<QVBoxLayout*> vLayouts;
    QVector<QPushButton*> pushButtons;
    QVector<QCPColorMap*> maps;
    QVector<QCPColorScale*> colorScales;
    QVector<QDoubleSpinBox*> SpinBoxes;
    QVector<QComboBox*> comboxes;

    // C-scan settings
    int CscanPlotMode=2; // default single slice

    // calculations
    void calculateSurface();

protected:
    // layout
    QVBoxLayout *layout=nullptr;

    // 3D dataset
    QVector<QVector<QVector<double>>> C_scan_double;
    QVector<QVector<QVector<std::complex<double>>>> C_scan_AS;

};

#endif // ULTRASOUND_CSCAN_PROCESS_H
