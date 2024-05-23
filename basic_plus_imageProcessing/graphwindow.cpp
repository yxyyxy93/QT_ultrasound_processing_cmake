#include "graphwindow.h"
#include <QApplication>
#include <QMainWindow>
#include <QVBoxLayout>
#include "qcustomplot.h"

GraphWindow::GraphWindow(QWidget *parent): QMainWindow(parent)
{
    centralWidget = new QWidget;
    plot = new QCustomPlot(centralWidget);
    colorMap = new QCPColorMap(plot->xAxis, plot->yAxis);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->addWidget(plot);
    centralWidget->setLayout(layout);
    setCentralWidget(centralWidget);
    resize(800, 600); // Adjust window size as needed

    // Hide the axes
    plot->xAxis->setVisible(false);
    plot->yAxis->setVisible(false);
    plot->xAxis2->setVisible(false); // Access the upper X axis, if used
    plot->yAxis2->setVisible(false); // Access the right Y axis, if used

    // // Create a color scale (color bar) for the color map
    // QCPColorScale *colorScale = new QCPColorScale(plot);
    // plot->plotLayout()->addElement(0, 1, colorScale); // Add color scale to the plot
    // colorMap->setColorScale(colorScale); // Tie the color map to the color scale
    // colorScale->setType(QCPAxis::atRight); // Position at right of plot

    // Set the Jet-like gradient
    QCPColorGradient gradient;
    gradient.loadPreset(QCPColorGradient::gpJet);
    colorMap->setGradient(gradient);

    plot->rescaleAxes();
    plot->replot();
}

GraphWindow::~GraphWindow()
{
    // Clean up dynamically allocated resources
}

void GraphWindow::setData(const QVector<QVector<double>>& data)
{
    if (!data.isEmpty() && !data.first().isEmpty())
    {
        int nx = data.size();
        int ny = data.first().size();
        colorMap->data()->setSize(nx, ny);  // Set the size of the color map data

        // Iterate over the data and set each cell
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                double value = data[i][j];
                if (value == -1) {
                    // Set transparent for -1 values
                    colorMap->data()->setCell(i, j, value);
                    colorMap->data()->setAlpha(i, j, 0);  // Make this cell fully transparent
                } else {
                    colorMap->data()->setCell(i, j, value);
                    colorMap->data()->setAlpha(i, j, 255);  // Fully opaque
                }
            }
        }

        colorMap->data()->setRange(QCPRange(0, nx), QCPRange(0, ny));
        plot->rescaleAxes();
        colorMap->rescaleDataRange();
        plot->replot();
    }
}


void GraphWindow::setColorScaleRange(double min, double max)
{
    colorMap->setDataRange(QCPRange(min, max));
    plot->replot();  // Redraw the plot with the new scale range
}

void GraphWindow::savePlotAsImage(const QString& filePath, const QString& format, int quality = -1)
{
    // Calculate image size based on the size of the data in the color map
    int nx = colorMap->data()->keySize();  // Number of columns in the data
    int ny = colorMap->data()->valueSize();  // Number of rows in the data

    // Determine a scaling factor to convert data dimensions to pixels
    int scaleFactor = 5;  // This is an arbitrary choice, adjust based on your needs

    // Calculate the image dimensions
    int imageWidth = nx * scaleFactor;
    int imageHeight = ny * scaleFactor;

    // Remove margins around the plot
    plot->axisRect()->setAutoMargins(QCP::msNone);
    plot->axisRect()->setMargins(QMargins(0, 0, 0, 0));

    if (format.toLower() == "png") {
        plot->savePng(filePath, imageWidth, imageHeight); // You can adjust the size as needed
    } else if (format.toLower() == "jpg" || format.toLower() == "jpeg") {
        plot->saveJpg(filePath, imageWidth, imageHeight, quality); // JPEG quality can be adjusted
    } else {
        qDebug() << "Unsupported file format";
    }
}
