// ImagePlotWidget.h

#ifndef IMAGEPLOTWIDGET_H
#define IMAGEPLOTWIDGET_H

#include <QWidget>
#include <QPainter>
#include <QVector>

class ImagePlotWidget : public QWidget {
    Q_OBJECT

public:
    explicit ImagePlotWidget(QWidget *parent = nullptr);
    void setData(const QVector<QVector<double>>& data);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    QVector<QVector<double>> data;
};

#endif // IMAGEPLOTWIDGET_H
