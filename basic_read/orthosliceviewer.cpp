#include "OrthosliceViewer.h"
#include "qcustomplot.h"
#include "utils.h"

#include <QDebug>

OrthosliceViewer::OrthosliceViewer(QWidget *parent,
                                   const QVector<QVector<QVector<std::complex<double>>>> &C_scan_AS,
                                   const QVector<QVector<QVector<double>>> &C_scan_double)
    : QWidget(parent), C_scan_AS(C_scan_AS), C_scan_double(C_scan_double), initialValueScrollBar(nullptr), decayRatioScrollBar(nullptr) {
    setupUI();
    connectSignals();
}

OrthosliceViewer::OrthosliceViewer(QWidget *parent,
                                   const QVector<QVector<QVector<double>>> &C_scan_double)
    : QWidget(parent), C_scan_double(C_scan_double), initialValueScrollBar(nullptr), decayRatioScrollBar(nullptr) {
    setupUI();
    connectSignals();
}

OrthosliceViewer::~OrthosliceViewer() {
    // Clean up dynamic memory, if any
}

void OrthosliceViewer::setdata(const QVector<QVector<QVector<std::complex<double>>>> &C_scan_AS,
                               const QVector<QVector<QVector<double>>> &C_scan_double) {
    this->C_scan_double = C_scan_double;
    this->C_scan_AS = C_scan_AS;
};

void OrthosliceViewer::setdata(const QVector<QVector<QVector<double>>> &C_scan_double) {
    this->C_scan_double = C_scan_double;
    this->C_scan_AS = C_scan_AS;
};

void OrthosliceViewer::setupUI() {
    // Create the QCustomPlot widgets
    customPlot1 = new QCustomPlot();
    customPlot2 = new QCustomPlot();
    customPlot3 = new QCustomPlot();
    customPlot_Ascan = new QCustomPlot(); // Initialize but do not add to layout yet
    // Create the ScrollBars and Labels
    scrollBarX = new QScrollBar(Qt::Horizontal);
    scrollBarY = new QScrollBar(Qt::Horizontal);
    scrollBarZ = new QScrollBar(Qt::Horizontal);
    sBX_label = new QLabel();
    sBY_label = new QLabel();
    sBZ_label = new QLabel();
    // Create the combo box and add the selection options
    comboBox = new QComboBox();
    comboBox->addItem("Max");
    comboBox->addItem("Average");
    comboBox->addItem("Single Index");
    comboBox->addItem("Single Index for original");
    // Create the main horizontal layout
    // Create the main horizontal layout
    QHBoxLayout *hLayout = new QHBoxLayout();
    qDebug() << "hLayout created";
    // Create the main vertical layout
    this->mainVLayout = new QVBoxLayout(this);
    qDebug() << "mainVLayout created";
    // Add the horizontal layout to the main vertical layout
    this->mainVLayout->addLayout(hLayout);
    this->mainVLayout->addWidget(customPlot_Ascan);
    // Create and populate vertical layouts for each plot
    QVBoxLayout *plot1l = new QVBoxLayout();
    plot1l->addWidget(this->scrollBarX);
    plot1l->addWidget(sBX_label);
    plot1l->addWidget(this->customPlot1);
    QVBoxLayout *plot2l = new QVBoxLayout();
    plot2l->addWidget(this->scrollBarY);
    plot2l->addWidget(sBY_label);
    plot2l->addWidget(this->customPlot2);
    QVBoxLayout *plot3l = new QVBoxLayout();
    // Add the vertical layouts to the main horizontal layout
    // Add the vertical layouts to the main horizontal layout
    hLayout->addLayout(plot1l);
    hLayout->addLayout(plot2l);
    hLayout->addLayout(plot3l);
    qDebug() << "Layouts added to hLayout";
    // Set the main layout for the widget
    this->setLayout(this->mainVLayout);
    qDebug() << "mainVLayout set as layout for the widget";
    // Set the range of the QScrollBars based on the size of the data
    this->scrollBarX->setRange(0, this->C_scan_AS.size() - 1);
    this->scrollBarY->setRange(0, this->C_scan_AS[0].size() - 1);
    this->scrollBarZ->setRange(0, this->C_scan_AS[0][0].size() - 1);
    // Set the initial values of the QScrollBars
    this->scrollBarX->setValue(this->C_scan_AS.size() / 2);
    this->scrollBarY->setValue(this->C_scan_AS[0].size() / 2);
    this->scrollBarZ->setValue(this->C_scan_AS[0][0].size() / 2);
    this->customPlot1->setMinimumSize(200, 200); // Example size, adjust as needed
    this->customPlot1->setMaximumHeight(1280); // Example size, adjust as needed
    this->customPlot2->setMinimumSize(200, 200); // Example size, adjust as needed
    this->customPlot2->setMaximumHeight(1280); // Example size, adjust as needed
    this->customPlot3->setMinimumSize(200, 200); // Example size, adjust as needed
    this->customPlot3->setMaximumHeight(1280); // Example size, adjust as needed
    // Modify vertical layouts
    plot1l->addWidget(scrollBarX);
    plot1l->addWidget(sBX_label);
    plot1l->addWidget(this->customPlot1);
    //
    plot2l->addWidget(scrollBarY);
    plot2l->addWidget(sBY_label);
    plot2l->addWidget(this->customPlot2);
    //
    plot3l->addWidget(scrollBarZ);
    plot3l->addWidget(sBZ_label);
    plot3l->addWidget(comboBox);
    plot3l->addWidget(this->customPlot3);
    // Initialize the gate and save button
    gateSaveButton = new QPushButton("Gate and Save Structure");
    this->mainVLayout->addWidget(gateSaveButton); // Add it to the main vertical layout
    connect(gateSaveButton, &QPushButton::clicked, this, &OrthosliceViewer::onGateSaveClicked);
    qDebug() << "finish add widgets";
    // After making changes to the layout, update the widget
    this->update(); // This will schedule a repaint event for the widget
    this->adjustSize(); // This will adjust the size of the widget to fit its contents
    this->connectSignals();
}

void OrthosliceViewer::connectSignals() {
    QObject::connect(this->scrollBarX, &QScrollBar::valueChanged, this,
                     &OrthosliceViewer::updatePlot);
    QObject::connect(this->scrollBarY, &QScrollBar::valueChanged, this,
                     &OrthosliceViewer::updatePlot);
    QObject::connect(this->scrollBarZ, &QScrollBar::valueChanged, this,
                     &OrthosliceViewer::updatePlot);
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarX, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarX->minimum()).arg(scrollBarX->maximum()).arg(value);
        sBX_label->setText(labelText);});
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarY, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarY->minimum()).arg(scrollBarY->maximum()).arg(value);
        sBY_label->setText(labelText);});
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarZ, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarZ->minimum()).arg(scrollBarZ->maximum()).arg(value);
        sBZ_label->setText(labelText);});

    // Connect the currentIndexChanged signal to a slot
    connect(comboBox, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &OrthosliceViewer::updateCscanPlotSelection);
}

void OrthosliceViewer::updatePlot() {
    int sliceX = this->scrollBarX->value();
    int sliceY = this->scrollBarY->value();
    int sliceZ = this->scrollBarZ->value();
    // Create QCPColorMap objects
    QCPColorMap *map1 = new QCPColorMap(this->customPlot1->xAxis,
                                        this->customPlot1->yAxis);
    QCPColorMap *map2 = new QCPColorMap(this->customPlot2->xAxis,
                                        this->customPlot2->yAxis);
    QCPColorMap *map3 = new QCPColorMap(this->customPlot3->xAxis,
                                        this->customPlot3->yAxis);
    // set the data for each QCPColorMap
    int x_size = this->C_scan_AS.size();
    int y_size = this->C_scan_AS[0].size();
    int z_size = this->C_scan_AS[0][0].size();
    map1->data()->setSize(y_size, z_size);
    map1->data()->setRange(QCPRange(0, y_size),
                           QCPRange(0, z_size));
    map2->data()->setSize(x_size, z_size);
    map2->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, z_size));
    map3->data()->setSize(x_size, y_size);
    map3->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    if (this->CscanPlotMode == 0) { // Max
        for (int i = 0; i < y_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map1->data()->setCell(i, j, std::abs(this->C_scan_AS[sliceX][i][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map2->data()->setCell(i, j, std::abs(this->C_scan_AS[i][sliceY][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < y_size; ++j) {
                double max_value = -std::numeric_limits<double>::max();
                for (int idx = std::max(sliceZ - 20, 0); idx < std::min(sliceZ + 20, z_size); ++idx)
                    max_value = std::max(max_value, std::abs(this->C_scan_AS[i][j][idx]));
                map3->data()->setCell(i, j, max_value);
            }
        }
    } else if (this->CscanPlotMode == 1) { // Average
        for (int i = 0; i < y_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map1->data()->setCell(i, j, std::abs(this->C_scan_AS[sliceX][i][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map2->data()->setCell(i, j, std::abs(this->C_scan_AS[i][sliceY][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < y_size; ++j) {
                double sum = 0;
                int count = 0;
                for (int idx = std::max(sliceZ - 5, 0); idx < std::min(sliceZ + 5, z_size); ++idx) {
                    sum += std::abs(this->C_scan_AS[i][j][idx]);
                    ++count;
                }
                double average = count > 0 ? sum / count : 0;
                map3->data()->setCell(i, j, average);
            }
        }
    } else if (this->CscanPlotMode == 2) { // Single Index
        for (int i = 0; i < y_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map1->data()->setCell(i, j, std::arg(this->C_scan_AS[sliceX][i][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map2->data()->setCell(i, j, std::arg(this->C_scan_AS[i][sliceY][j]));
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < y_size; ++j) {
                map3->data()->setCell(i, j, std::arg(this->C_scan_AS[i][j][sliceZ]));
            }
        }
    } else if (this->CscanPlotMode == 3) { // Original
        for (int i = 0; i < y_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map1->data()->setCell(i, j, this->C_scan_double[sliceX][i][j]);
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < z_size; ++j) {
                map2->data()->setCell(i, j, this->C_scan_double[i][sliceY][j]);
            }
        }
        for (int i = 0; i < x_size; ++i) {
            for (int j = 0; j < y_size; ++j) {
                map3->data()->setCell(i, j, this->C_scan_double[i][j][sliceZ]);
            }
        }
    }
    map3->setGradient(QCPColorGradient::gpJet);
    // Rescale the color map data range to fit the new data
    map1->rescaleDataRange(true);
    map2->rescaleDataRange(true);
    map3->rescaleDataRange(true);
    map1->rescaleAxes();
    map2->rescaleAxes();
    map3->rescaleAxes();
    this->customPlot1->xAxis->setRange(0, y_size);
    this->customPlot1->yAxis->setRange(0, z_size);
    this->customPlot2->xAxis->setRange(0, x_size);
    this->customPlot2->yAxis->setRange(0, z_size);
    this->customPlot3->xAxis->setRange(0, x_size);
    this->customPlot3->yAxis->setRange(0, y_size);
    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    this->customPlot3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    // Call replot() to update the plot with the new data
    this->customPlot1->replot();
    this->customPlot2->replot();
    this->customPlot3->replot();
    // Connect the plot's mouse press signal to the onCustomPlotClicked slot
    connect(this->customPlot3, &QCustomPlot::mousePress,
            this, &OrthosliceViewer::onCustomPlotClicked_Cscan);
}

void OrthosliceViewer::onCustomPlotClicked_Cscan(QMouseEvent* event) {
    static QElapsedTimer lastClickTime;
    int minClickInterval = 500; // minimum time between clicks in milliseconds
    //
    if (event->button() == Qt::LeftButton)
    {
        if (lastClickTime.isValid() && lastClickTime.elapsed() < minClickInterval)
        {
            // ignore click if it's too soon after the last one
            return;
        }
        lastClickTime.start();
        // handle the click event
        int x = this->customPlot3->xAxis->pixelToCoord(event->pos().x());
        int y = this->customPlot3->yAxis->pixelToCoord(event->pos().y());
        qDebug() << "Clicked at (" << x << "," << y << ")";
        // plot Ascan
        this->signal = this->C_scan_double[x][y];
        for (int i = 0; i < this->signal.size(); ++i) {
            this->time.append(i);
        }
        // calculate the analytic-signal
        QVector<std::complex<double>> Ascan_as;
        Ascan_as = analyticSignal(signal);
        QVector<double> Ascan_as_abs;
        for (const auto& element : Ascan_as) {
            double absoluteValue = std::abs(element);
            Ascan_as_abs.append(absoluteValue);
        }
        // create a new
        if (this->customPlot_Ascan == nullptr) {
            this->customPlot_Ascan = new QCustomPlot(this);
            QVBoxLayout *mainLayout = dynamic_cast<QVBoxLayout*>(this->layout());
            if (mainLayout) {
                mainLayout->addWidget(this->customPlot_Ascan);
            }
        }
        this->customPlot_Ascan->clearGraphs();
        // add two new graphs and set their look:
        this->customPlot_Ascan->addGraph();
        this->customPlot_Ascan->graph(0)->setPen(QPen(Qt::blue)); // line color blue for first graph
        this->customPlot_Ascan->graph(0)->setBrush(QBrush(QColor(0, 0, 255, 20))); // first graph will be filled with translucent blue
        this->customPlot_Ascan->addGraph();
        this->customPlot_Ascan->graph(1)->setPen(QPen(Qt::red)); // line color red for second graph
        // configure right and top axis to show ticks but no labels:
        // (see QCPAxisRect::setupFullAxesBox for a quicker method to do this)
        this->customPlot_Ascan->xAxis2->setVisible(true);
        this->customPlot_Ascan->xAxis2->setTickLabels(false);
        this->customPlot_Ascan->yAxis2->setVisible(true);
        this->customPlot_Ascan->yAxis2->setTickLabels(false);
        // make left and bottom axes always transfer their ranges to right and top axes:
        connect(this->customPlot_Ascan->xAxis,
                SIGNAL(rangeChanged(QCPRange)),
                this->customPlot_Ascan->xAxis2,
                SLOT(setRange(QCPRange)));
        connect(this->customPlot_Ascan->yAxis,
                SIGNAL(rangeChanged(QCPRange)),
                this->customPlot_Ascan->yAxis2,
                SLOT(setRange(QCPRange)));
        // pass data points to graphs:
        this->customPlot_Ascan->graph(0)->setData(this->time,
                                                  this->signal);
        // this->customPlot_Ascan->graph(1)->setData(time,
        //                                           Ascan_as_abs);
        // Generate the Exponential Decay Vector
        int decayVectorSize = this->signal.size(); // Matching the size of the signal vector
        // Find the maximum of Ascan_as_abs
        double maxValue = *std::max_element(Ascan_as_abs.begin(), Ascan_as_abs.end());
        double decayRate = 0.5; // Example decay rate
        decayVector.clear();
        decayVector.resize(this->signal.size()); // Set the size of decayVector to match meanAscan
        // Initialize vectors to store mean and std values, assuming all ascans have the same size
        QVector<double> meanAscan(this->C_scan_double.first().first().size(), 0);
        QVector<double> stdAscan(this->C_scan_double.first().first().size(), 0);
        // Calculate the sum for mean calculation
        for (int timeInstant = 0; timeInstant < meanAscan.size(); ++timeInstant) {
            double sum = 0;
            int count = 0;
            for (const auto& xyPlane : this->C_scan_double) {
                for (const auto& xLine : xyPlane) {
                    sum += abs(xLine[timeInstant]); // Ensure abs() is valid for your data type
                    ++count;
                }
            }
            meanAscan[timeInstant] = sum / count;
        }
        // Step 2: Calculate Standard Deviation for Each Time Instant
        for (int timeInstant = 0; timeInstant < meanAscan.size(); ++timeInstant) {
            double varianceSum = 0;
            int count = 0;
            for (const auto& xyPlane : this->C_scan_double) {
                for (const auto& xLine : xyPlane) {
                    double diff = abs(xLine[timeInstant]) - meanAscan[timeInstant];
                    varianceSum += diff * diff;
                    ++count;
                }
            }
            double variance = varianceSum / count;
            stdAscan[timeInstant] = sqrt(variance);
            decayVector[timeInstant] = meanAscan[timeInstant] + 1 * stdAscan[timeInstant];
        }
        // Add a third graph for the Exponential Decay Vector
        this->customPlot_Ascan->addGraph();
        // Create a QPen object with the color green
        QPen pen(Qt::green);
        // Set the pen width to a larger size, for example, 3 pixels
        pen.setWidth(5);
        // Apply the pen to the third graph
        this->customPlot_Ascan->graph(2)->setPen(pen);
        this->customPlot_Ascan->graph(2)->setData(this->time, decayVector); // Use the same 'time' vector for x-axis
        this->customPlot_Ascan->rescaleAxes();
        // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
        this->customPlot_Ascan->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        this->customPlot_Ascan->replot();
        // Assuming 'this' is a QWidget or QMainWindow derivative
        if (!initialValueScrollBar) {
            initialValueScrollBar = new QScrollBar(Qt::Horizontal, this);
            initialValueScrollBar->setRange(1, 100); // Adjust range as needed
            connect(initialValueScrollBar, &QScrollBar::valueChanged, this, [this, maxValue](int newValue) {
                this->onInitialValueChanged(newValue, maxValue);
            });
        }
        if (!decayRatioScrollBar) {
            decayRatioScrollBar = new QScrollBar(Qt::Horizontal, this);
            decayRatioScrollBar->setRange(1, 100); // Decay ratio in percentage (1-100)
            connect(decayRatioScrollBar, &QScrollBar::valueChanged, this, [this, maxValue](int newValue) {
                this->onDecayRatioChanged(newValue, maxValue);
            });
        }
        // Ensure scrollbars are not already part of any layout
        bool scrollBarsAlreadyAdded = false;
        QLayout *parentLayout = initialValueScrollBar->parentWidget()->layout();
        if (parentLayout) {
            for (int i = 0; i < parentLayout->count(); ++i) {
                QWidget *widget = parentLayout->itemAt(i)->widget();
                if (widget == initialValueScrollBar || widget == decayRatioScrollBar) {
                    scrollBarsAlreadyAdded = true;
                    break;
                }
            }
        }
        if (!scrollBarsAlreadyAdded) {
            QHBoxLayout *hScrollBarLayout = new QHBoxLayout();
            hScrollBarLayout->addWidget(initialValueScrollBar);
            hScrollBarLayout->addWidget(decayRatioScrollBar);
            this->mainVLayout->addLayout(hScrollBarLayout);
        }
    }
}

// *********** for C-scan settings
void OrthosliceViewer::updateCscanPlotSelection(int index) {
    this->CscanPlotMode = index; // This should be a member variable to store the current selection
    qDebug() << index;
}

// ************* for time gate
void OrthosliceViewer::onInitialValueChanged(int value, double max_value) {
    currentInitialValue = static_cast<double>(value)/50;
    // Regenerate the exponential decay vector
    QVector<double> decayVector = createExponentialDecayVector(this->signal.size(), currentInitialValue * max_value, currentDecayRate);
    // Update the graph
    this->customPlot_Ascan->graph(2)->setData(this->time, decayVector);
    this->customPlot_Ascan->replot();
}

void OrthosliceViewer::onDecayRatioChanged(int value, double max_value) {
    currentDecayRate = static_cast<double>(value) / 100.0; // Convert to a fraction
    // Regenerate the exponential decay vector
    QVector<double> decayVector = createExponentialDecayVector(this->signal.size(), max_value, currentDecayRate);
    // Update the graph
    this->customPlot_Ascan->graph(2)->setData(time, decayVector);
    this->customPlot_Ascan->replot();
}

void OrthosliceViewer::onGateSaveClicked() {
    // Use QFileDialog to get the file name from the user
    QString fileName = QFileDialog::getSaveFileName(this, tr("Save Data"), "", tr("CSV Files (*.csv);;All Files (*)"));
    // Check if the user canceled the dialog
    if (fileName.isEmpty()) {
        return; // Exit if no file name was provided
    }

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        qDebug() << "Error opening file for writing";
        return;
    }

    QTextStream out(&file);
    out << this->C_scan_double.size() << "," << this->C_scan_double[0].size() << "," << this->C_scan_double[0][0].size() << "\n";
    for (int i = 0; i < this->C_scan_double.size(); ++i) {
        for (int j = 0; j < this->C_scan_double[i].size(); ++j) { // Ensuring to use dynamic size of second dimension
            const QVector<double>& Ascan = this->C_scan_double[i][j];
            QStringList values;
            for (int k = 0; k < Ascan.size(); ++k) {
                double Ascan_abs = std::abs(Ascan[k]); // Correctly computing the absolute value
                if (Ascan_abs < 2*this->decayVector[k]) { // Corrected indexing
                    values << QString::number(0);
                } else {
                    // values << QString::number((Ascan_abs - this->decayVector[k])/Ascan_abs);
                    values << QString::number(20 * log10(Ascan_abs/this->decayVector[k]));
                }
            }
            out << values.join(',') << "\n";
        }
    }
    file.close();
    qDebug() << "Data saved to gatedData.csv";
}

