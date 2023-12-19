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
    resize(800, 600);
}

GraphWindow::~GraphWindow()
{

}

void GraphWindow::setData(QVector<QVector<double>> data)
{
    // Check if data is not empty and has at least one vector
    if (!data.empty() && !data[0].empty())
    {
        int nx = data.size(); // number of columns
        int ny = data[0].size(); // number of rows

        colorMap->data()->setSize(nx, ny);

        // Fill out the data
        for (int i = 0; i < nx; ++i)
            for (int j = 0; j < ny; ++j)
                colorMap->data()->setCell(i, j, data[i][j]);

        // Set the range according to the data dimension
        colorMap->data()->setRange(QCPRange(0, nx), QCPRange(0, ny));

        // Rescale data dimension and data range to fit
        plot->rescaleAxes();
        colorMap->rescaleDataRange();
    }
}
