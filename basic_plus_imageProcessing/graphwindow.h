#ifndef GRAPHWINDOW_H
#define GRAPHWINDOW_H

#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include "qcustomplot.h"

class GraphWindow : public QMainWindow
{
    Q_OBJECT
public:
    GraphWindow(QWidget *parent = nullptr);
    ~GraphWindow();
    void setData(QVector<QVector<double>> data);
private:
    QWidget *centralWidget;
    QCustomPlot *plot;
    QCPColorMap *colorMap;
};

#endif // GRAPHWINDOW_H
