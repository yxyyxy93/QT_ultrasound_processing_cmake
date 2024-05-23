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
    void setData(const QVector<QVector<double>> &data);
    void setColorScaleRange(double min, double max);
    void savePlotAsImage(const QString& filePath, const QString& format, int quality);

private:
    QWidget *centralWidget;
    QCustomPlot *plot;
    QCPColorMap *colorMap;
};

#endif // GRAPHWINDOW_H
