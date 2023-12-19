// ImagePlotWidget.cpp

#include "ImagePlotWidget.h"

ImagePlotWidget::ImagePlotWidget(QWidget *parent) : QWidget(parent) {}

void ImagePlotWidget::setData(const QVector<QVector<double>>& data) {
    this->data = data;
    this->update();  // Trigger a repaint
}

void ImagePlotWidget::paintEvent(QPaintEvent *event) {
    QPainter painter(this);
    if (data.isEmpty()) return;

    int rows = data.size();
    int cols = data[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Map the double value to a grayscale value (0-255)
            int grayScale = static_cast<int>(data[i][j] * 255);
            grayScale = qBound(0, grayScale, 255); // Ensure value is within 0-255

            QColor color(grayScale, grayScale, grayScale);
            painter.setPen(color);
            painter.drawPoint(j, i); // Plot each point
        }
    }
}
