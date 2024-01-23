#ifndef ORTHOSLICEVIEWER_H
#define ORTHOSLICEVIEWER_H

#include <QWidget>
#include <QCustomPlot.h>
#include <QScrollBar>
#include <QLabel>
#include <QComboBox>
#include <QPushButton>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <complex>

class OrthosliceViewer : public QWidget {
    Q_OBJECT

public:
    explicit OrthosliceViewer(QWidget *parent = nullptr,
                              const QVector<QVector<QVector<std::complex<double>>>> &C_scan_AS = {},
                              const QVector<QVector<QVector<double>>> &C_scan_double = {});

    explicit OrthosliceViewer(QWidget *parent = nullptr,
                              const QVector<QVector<QVector<double>>> &C_scan_double = {});

    ~OrthosliceViewer();

    void setdata(const QVector<QVector<QVector<std::complex<double>>>> &C_scan_AS,
                const QVector<QVector<QVector<double>>> &C_scan_double);
    void setdata(const QVector<QVector<QVector<double>>> &C_scan_double);

    void updatePlot();
    void onCustomPlotClicked_Cscan(QMouseEvent* event);

    void updateCscanPlotSelection(int index);

private:
    QCustomPlot *customPlot1;
    QCustomPlot *customPlot2;
    QCustomPlot *customPlot3;
    QCustomPlot *customPlot_Ascan;
    QScrollBar *scrollBarX;
    QScrollBar *scrollBarY;
    QScrollBar *scrollBarZ;
    QLabel *sBX_label;
    QLabel *sBY_label;
    QLabel *sBZ_label;
    QComboBox *comboBox;
    QHBoxLayout *hLayout;
    QVBoxLayout *plot1l;
    QVBoxLayout *plot2l;
    QVBoxLayout *plot3l;

    QVector<QVector<QVector<std::complex<double>>>> C_scan_AS;
    QVector<QVector<QVector<double>>> C_scan_double;
    int CscanPlotMode;

    void setupUI();
    void connectSignals();

};

#endif // ORTHOSLICEVIEWER_H
