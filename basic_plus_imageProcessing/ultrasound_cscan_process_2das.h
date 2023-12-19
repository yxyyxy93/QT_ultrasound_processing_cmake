#ifndef ULTRASOUND_CSCAN_SEG_H
#define ULTRASOUND_CSCAN_SEG_H

#include <QMainWindow>
#include <QTimer>
#include <QtWidgets>
#include "../basic_read/ultrasound_cscan_process.h"
#include <memory>
#include "imageplotwidget.h"

// do not need to set default arguments for a constructor in a derived class
// if the base class already has default arguments set in its constructor.
class ultrasound_cscan_process_2DAS: public ultrasound_Cscan_process
{

public:
    ultrasound_cscan_process_2DAS(QWidget *parent,
                                  QString fn,
                                  int fs,
                                  double fx,
                                  double fy);
    ~ultrasound_cscan_process_2DAS();

private:
    QVBoxLayout *layout2;

    QProgressBar *m_progressBar_page2;

    // Qt plots
    QCustomPlot *customPlot1_page2;
    QCustomPlot *customPlot2_page2;
    QCustomPlot *customPlot3_page2;
    QScrollBar *scrollBarZ_page2;

    QLabel* imageLabel;

    // 2D AS properties
    QVector<QVector<double>> img_Cscan;
    double sc;
    double sf;
    int n;

     // ******** 3rd page
    QWidget *page3;
    QVBoxLayout *layout3;
    QCustomPlot *customPlot_Cscan;
    QCustomPlot *customPlot_structure;
    int lastValueX;
    int lastValueY;
    QTimer *debounceTimer_Cscan;
    QTimer *debounceTimer_structure;

    QVector<QVector<QVector<double>>> structure;

public slots:
    // ******** 3rd page
    void onCscanScrollBarChanged(int value);
    void onStructureScrollBarChanged(int value);
    void updateCscanPlot();
    void updateStructurePlot();
    void SavePlot();
    void onPlot3DButtonClicked();

private slots:
    // 2D　analytic-signal
    void handleButton_2DAS();
    void updatePlot_page2();
    void handleButton_2DAS_oneslice();

    // in "handleButton_surface"
    void handleButton_myButton_showkernel();

    // statistic for 2d as
    void handleButton_myButton_stat();

    // ******** 3rd page
    void load_structures();
    void Plot_Ccans_structures();
};

#endif // ULTRASOUND_CSCAN_PROCESS_2DAS_H
