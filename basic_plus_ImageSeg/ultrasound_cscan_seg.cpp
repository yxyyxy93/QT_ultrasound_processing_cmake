#include "ultrasound_cscan_seg.h"
//#include "ui_ultrasound_cscan_seg.h"

#include <QVBoxLayout>
#include <QMenuBar>
#include <QMenu>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QTextEdit>
#include <customgraphicsscene.h>
#include <complex>
#include "JetColorMap.h"
#include "..\basic_read\utils.h"

ultrasound_cscan_seg::ultrasound_cscan_seg(QWidget *parent,
                                           QString fn,
                                           int fs,
                                           double fx,
                                           double fy)
    : ultrasound_Cscan_process(parent, fn, fs, fx, fy)
{
    // Create the tab widget
    QTabWidget *tabWidget = new QTabWidget();
    // ********************* First page
    tabWidget->addTab(this, "Simple C-scan check 1");
    // ******************** 2nd page
    QWidget *page2 = new QWidget();
    this->layout2 = new QVBoxLayout(page2);
    QPushButton *myButton_Cscan = new QPushButton(tr("Plot C-scan images"), page2);
    connect(myButton_Cscan, &QPushButton::clicked, this, &ultrasound_cscan_seg::handleButton_Cscan);
    this->layout2->addWidget(myButton_Cscan);
    tabWidget->addTab(page2, "Image segmentation");
    // ************* 3rd page
    QWidget *page3 = new QWidget();
    this->layout3 = new QVBoxLayout(page3);
    tabWidget->addTab(page3, "Noise add and save");
    // Show the tab widget
    tabWidget->show();
    // Initialize the progress bar for page3
    this->progressBarPage3 = new QProgressBar(page3);
    this->progressBarPage3->setRange(0, 100);
    this->layout3->addWidget(this->progressBarPage3);
    // select folder
    QPushButton *selectFolderButton = new QPushButton(tr("select folder"), page3);
    this->folderLabel = new QLabel(page3);
    // Create a horizontal layout
    QHBoxLayout *rowLayout = new QHBoxLayout(page3);
    // Add the button and label to the horizontal layout
    rowLayout->addWidget(selectFolderButton);
    rowLayout->addWidget(this->folderLabel);
    rowLayout->setSpacing(10); // Adjust spacing as needed
    // Add the horizontal layout to your existing vertical layout (layout3)
    this->layout3->addLayout(rowLayout);
    connect(selectFolderButton, &QPushButton::clicked, this, &ultrasound_cscan_seg::selectFolder);
    // ************
    // Check if the 'Add Noise' button already exists in the layout
    if (!widgetExistsInLayout<QPushButton>(this->layout3, "process(add noise, downsample, corp) and save")) {
        // Create the 'Add Noise' button
        QPushButton* myButton_multiSNR = new QPushButton(tr("Add multi SNR (separated by ,)"), this);
        myButton_multiSNR->setObjectName("myButton_multiplySNR"); // Assign unique object name
        this->layout3->addWidget(myButton_multiSNR);
        connect(myButton_multiSNR, &QPushButton::clicked, this, &ultrasound_cscan_seg::handleButton_multiSNR);
    }
    // Check if the SNR input field already exists in the layout
    if (!widgetExistsInLayout<QLineEdit>(this->layout3, "multisnrInput")) {
        // Create the SNR input field
        // Create the SNR input field
        this->multisnrInput = new QLineEdit(this);
        this->multisnrInput->setObjectName("multisnrInput"); // Assign unique object name
        this->multisnrInput->setPlaceholderText("Enter multi SNR values");
        this->layout3->addWidget(this->multisnrInput);
    }
    // add a button to downsapmle the data
    if (!widgetExistsInLayout<QLineEdit>(this->layout3, "downsampleRateInput")) {
        // Create the Downsample Rate Input Field
        this->downsampleRateInput = new QLineEdit(this);
        this->downsampleRateInput->setObjectName("downsampleRateInput");
        this->downsampleRateInput->setPlaceholderText("Enter downsample rate");
        this->layout3->addWidget(this->downsampleRateInput);
    }
    // Check if the Crop Signal Input Field already exists in the layout
    if (!widgetExistsInLayout<QLineEdit>(this->layout3, "cropSignalInput")) {
        // Create the Crop Signal Input Field
        this->cropSignalInput = new QLineEdit(this);
        this->cropSignalInput->setObjectName("cropSignalInput"); // Assign unique object name
        this->cropSignalInput->setPlaceholderText("Enter crop values (start, end)");
        this->layout3->addWidget(this->cropSignalInput);
    }
}

// ************* 2D analytic-signal and visualization
void ultrasound_cscan_seg::handleButton_Cscan() {
    this->x_size = ultrasound_Cscan_process::C_scan_AS.size();
    this->y_size = (x_size>0) ? ultrasound_Cscan_process::C_scan_AS[0].size() : 0;
    this->z_size = (y_size>0) ? ultrasound_Cscan_process::C_scan_AS[0][0].size() : 0;
    // initialize c_scan_mask
    C_scan_mask = QVector<QVector<QVector<bool>>>(x_size, QVector<QVector<bool>>(y_size, QVector<bool>(z_size, false)));
    qDebug() << x_size << y_size << z_size;

    //    // ****************** create the orthoslice visual
    // Qcustomplot
    this->scrollBarZ_page2 = new QScrollBar(Qt::Horizontal);
    this->scrollBarZ_page2->setGeometry(50, 470, 600, 20);
    // add an value label
    QLabel *sBZ_label = new QLabel();
    sBZ_label->setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
    //    Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarZ_page2, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarZ_page2->minimum()).arg(scrollBarZ_page2->maximum()).arg(value);
        sBZ_label->setText(labelText);
    });
    //    // Create the QCustomPlot widget
    //    this->customPlot1_page2 = new QCustomPlot();
    //    // Create a horizontal layout for the main window
    QWidget *rightPanel = new QWidget();
    QHBoxLayout *hLayout = new QHBoxLayout();
    rightPanel->setLayout(hLayout);
    //    // Create a vertical layout for 1
    QVBoxLayout *plot1l = new QVBoxLayout();
    QWidget *plot1w = new QWidget();
    plot1w->setLayout(plot1l);
    //    // Add some widgets to the right panel
    //    plot1l->addWidget(this->customPlot1_page2);
    plot1l->addWidget(this->scrollBarZ_page2);
    plot1l->addWidget(sBZ_label);
    //    // Create a vertical layout for 2
    //    QVBoxLayout *plot2l = new QVBoxLayout();
    //    QWidget *plot2w = new QWidget();
    //    plot2w->setLayout(plot2l);
    //    // Add some widgets to the right panel
    //    plot2l->addWidget(this->customPlot2_page2);
    //    // Add some widgets to the right panel;
    hLayout->addWidget(plot1w);
    //    hLayout->addWidget(plot2w);
    this->layout2->addWidget(rightPanel);

    // Set the range of the QScrollBars based on the size of the data
    this->scrollBarZ_page2->setRange(0, z_size - 1);
    // Set the initial values of the QScrollBars
    this->scrollBarZ_page2->setValue(z_size / 2);
    // Connect the valueChanged() signals of the QScrollBars to update the plot data
    QTimer *debounceTimer = new QTimer(this); // avoid crash. won't get called until a short period after the user stops dragging the scrollbar.
    debounceTimer->setSingleShot(true);
    QObject::connect(this->scrollBarZ_page2, &QScrollBar::valueChanged, [this, debounceTimer]() {
        debounceTimer->start(100);  // 100 ms debounce
    });
    //    Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    QObject::connect(debounceTimer, &QTimer::timeout, [this, sBZ_label]() {
        int value = this->scrollBarZ_page2->value();
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarZ_page2->minimum()).arg(scrollBarZ_page2->maximum()).arg(value);
        sBZ_label->setText(labelText);
        //        this->updatePlot_page2();
    });
    Q_ASSERT(this->scrollBarZ_page2 != nullptr);

    // ************** push botton to draw lines
    QPushButton *myButton_draw = new QPushButton(tr("Plot and draw lines"), this);
    connect(myButton_draw, &QPushButton::clicked, this,
            &ultrasound_cscan_seg::handleButton_draw);
    this->layout2->addWidget(myButton_draw);
    // ************** add a embeded window
    this->embeddedView = new QGraphicsView(this);
    this->layout2->addWidget(embeddedView);
}

// ********************** on page 2 segmentation manipulation
void ultrasound_cscan_seg::handleButton_draw(){
    int sliceZ = this->scrollBarZ_page2->value();
    QImage image = convertToImage<std::complex<double>>(ultrasound_Cscan_process::C_scan_AS, sliceZ);

    if (this->scene == nullptr) {
        this->scene = new CustomGraphicsScene(this);
    }

    // Clear drawn lines in the CustomGraphicsScene
    this->scene->clearDrawnLines();
    this->scene->removeAllPolygons();

    //    connect(this->scene, &CustomGraphicsScene::lineDrawn, this, &ultrasound_cscan_seg::onLineDrawn);
    // Set the scene for the embedded view
    this->embeddedView->setScene(scene);
    //    scene->setSceneRect(0, 0, 800, 600);  // adjust the size as per your requirement
    scene->addPixmap(QPixmap::fromImage(image));
    this->embeddedView->show();

    // add buttons to draw masks and save multitudes of cscans
    // Check and add "Close Area" button
    if (!widgetExistsInLayout<QPushButton>(this->layout2, "Close Area")) {
        QPushButton *closeAreaButton = new QPushButton("Close Area", this);
        this->layout2->addWidget(closeAreaButton);
        connect(closeAreaButton, &QPushButton::clicked, this, &ultrasound_cscan_seg::closeDrawnArea);
    }
    // Check and add "Save multitude of images" button
    if (!widgetExistsInLayout<QPushButton>(this->layout2, "Save multitude of images")) {
        QPushButton *saveButton = new QPushButton("Save multitude of images", this);
        this->layout2->addWidget(saveButton);
        connect(saveButton, &QPushButton::clicked, this, &ultrasound_cscan_seg::saveImage);
    }
};

void ultrasound_cscan_seg::onLineDrawn(const QLineF &line) {
    // Handle the line as you wish, e.g., store it, process it, etc.
}

void ultrasound_cscan_seg::saveImage() {
    // Use QInputDialog to get the range values from the user
    bool ok1, ok2;
    int startValue = QInputDialog::getInt(this, tr("Input Start Value"),
                                          tr("Start Z-Index:"), 300, 0, 10000, 1, &ok1);
    int endValue = QInputDialog::getInt(this, tr("Input End Value"),
                                        tr("End Z-Index:"), 900, 0, 10000, 1, &ok2);
    int step = QInputDialog::getInt(this, tr("Input step Value"),
                                    tr("Step Z-Index:"), 10, 0, 1000, 1, &ok2);

    // If the user pressed Cancel, then ok1 or ok2 will be false
    if (!ok1 || !ok2) {
        qDebug() << "Input dialog canceled";
        return;
    }

    QString folderPath = QFileDialog::getExistingDirectory(this, tr("Select Root Directory to Save Images"));

    if (!folderPath.isEmpty()) {
        // Create subfolders
        QDir rootDir(folderPath);
        rootDir.mkdir("original_images");
        rootDir.mkdir("groundtruth");

        // Specify the paths to the subfolders
        QString subfolder1Path = folderPath + "/original_images";
        QString subfolder2Path = folderPath + "/groundtruth";

        // Assuming you have two QVector<QVector<QVector<double>>> members named 'data1' and 'data2'
        for(int z = startValue; z < endValue; z+=step) {
            QImage image1 = convertToImage<std::complex<double>>(ultrasound_Cscan_process::C_scan_AS, z);
            QImage image2 = convertToImage<bool>(this->C_scan_mask, z);

            QString image1FullPath = subfolder1Path + "/cscan_" + QString::number(z) + ".png";
            QString image2FullPath = subfolder2Path + "/cscan_mask" + QString::number(z) + ".png";

            image1.save(image1FullPath);
            image2.save(image2FullPath);
        }
    }

    qDebug() << "images saved";
}

void ultrasound_cscan_seg::closeDrawnArea() {
    const QVector<QLineF>& drawnLines = static_cast<CustomGraphicsScene*>(this->embeddedView->scene())->getDrawnLines();

    if (drawnLines.size() > 1) {
        // Close the area by connecting the last point to the first point
        QLineF closingLine(drawnLines.last().p2(), drawnLines.first().p1());

        // Convert the lines to a polygon and add it to the scene
        QPolygonF polygon;
        for (const QLineF &line : drawnLines) {
            polygon << line.p1();
        }
        polygon << closingLine.p1() << closingLine.p2();  // Add the closing line's points

        CustomGraphicsScene* scene = static_cast<CustomGraphicsScene*>(this->embeddedView->scene());
        scene->addPolygon(polygon, QPen(Qt::blue), QBrush(Qt::lightGray));

        // Create a mask image
        QRectF rect = this->embeddedView->scene()->sceneRect();
        QImage maskImage(rect.width(), rect.height(), QImage::Format_Grayscale8);
        maskImage.fill(Qt::black);
        QPainter painter(&maskImage);
        painter.setBrush(Qt::white);
        painter.setPen(Qt::white);
        painter.drawPolygon(polygon);
        painter.end();
        // Save the mask image
        maskImage.save("mask.png");

        // save to mask data
        int sliceZ = this->scrollBarZ_page2->value();
        for (int x = 0; x < x_size; x++) {
            for (int y = 0; y < y_size; y++) {
                // Check if the current point is inside the polygon
                for (int z = z_size-1; z>=sliceZ; z--)
                    if (polygon.containsPoint(QPointF(x, y), Qt::OddEvenFill)) {
                        this->C_scan_mask[x][y][z] = true;
                    }
            }
        }

        qDebug() << "mask data saved";
    }
}

template <typename T>
QImage ultrasound_cscan_seg::convertToImage(const QVector<QVector<QVector<T>>>& data, int z) {
    int x_size = data.size();
    int y_size = data[0].size();
    QImage image(x_size, y_size, QImage::Format_RGB32); // Adjust format as necessary

    // ********** find the maximum and minimum
    double maxValue = std::abs(data[0][0][z]);  // Initialize with first element
    double minValue = std::abs(data[0][0][z]);  // Initialize with first element
    for(const QVector<QVector<T>>& outerVector : data) {
        for (const QVector<T>& innerVector: outerVector){
            T value = innerVector[z];
            if(std::abs(value) > maxValue)
                maxValue = std::abs(value);
            if(std::abs(value) < minValue)
                minValue = std::abs(value);
        }
    }

    JetColorMap colormap;
    for(int x = 0; x < x_size; ++x) {
        for(int y = 0; y < y_size; ++y) {
            double value = std::abs(data[x][y][z]);
            // Adjust the normalization/translation of data value to pixel value as necessary
            double normalizedValue = (value - minValue) / (maxValue - minValue);
            QColor jetColor = colormap.valueToJet(normalizedValue);
            image.setPixelColor(x, y, jetColor);
        }
    }
    return image;
}

ultrasound_cscan_seg::~ultrasound_cscan_seg()
{
    //    delete ui;
}

// ****************** on page 3 ***********
void ultrasound_cscan_seg::handleButton_multiSNR() {
    QLineEdit* SNRInput = this->multisnrInput;
    qDebug() << SNRInput->text();
    // read downsample ratio
    bool ok_ds;
    int downsampleRate = this->downsampleRateInput->text().toInt(&ok_ds);
    if (!ok_ds || downsampleRate <= 0) {
        // Handle error: invalid input
        return;
    }
    // read the corping values for signal
    int startValue, endValue;
    this->readCropValues(startValue, endValue);
    // Check if values are valid, assuming that startValue should be less than endValue
    if (startValue < 0 || endValue <= startValue) {
        qDebug() << "Invalid cropping range. Start should be less than End.";
        return;
    }
    if (SNRInput) {
        QStringList snrValuesStr = SNRInput->text().split(',');
        foreach (const QString &snrStr, snrValuesStr) {  // Process each SNR value
            bool ok;
            double snrDb = snrStr.trimmed().toDouble(&ok);
            if (ok) {
                QString filename = this->fn;
                // Remove '.csv' from the base name
                int dotIndex = filename.indexOf('.'); // Find the index of the first dot
                if (dotIndex != -1) { // Check if a dot was found
                    filename = filename.left(dotIndex);
                    // Use modifiedString as needed
                } else {
                    // Handle the case where there's no dot in the string
                }
                // Create the directory
                QDir dir;
                if (!dir.mkpath(filename)) {
                    qWarning() << "Could not create directory:" << filename;
                } else {
                    qDebug() << "Directory created:" << filename;
                }
                // Create the new filename by appending the double value
                QString newFileName = filename + "/_snr_" + QString::number(snrDb, 'f', 2) + ".csv"; // 'f': fixed-point notation, '2': two decimal places
                qDebug() << "New file path:" << newFileName;
                QFile file(newFileName);
                if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                    qWarning() << "Could not open file for writing:" << file.errorString();
                    return;
                }
                QTextStream out(&file);
                // Write dimensions as the first line
                if (!this->C_scan_double.isEmpty() && !this->C_scan_double[0].isEmpty()) {
                    int x = this->C_scan_double.size();
                    int y = this->C_scan_double[0].size();
                    int z = endValue-startValue+1;
                    out << x << "," << y << "," << z << "\n";
                } else {
                    qWarning() << "Data is empty, cannot write to file";
                    file.close();
                    return;
                }
                for (int i = 0; i < this->C_scan_double.size(); i++) {
                    for (int j = 0; j < this->C_scan_double[i].size(); j++) {
                        QVector<double> Ascan = this->C_scan_double[i][j]; // avoid too small value
                        addGaussianNoise(Ascan, snrDb);
                        QVector<double> Ascan_ds = downsampleVector(Ascan, downsampleRate);
                        QStringList values;
                        for (int k = startValue; k <= endValue; ++k) {
                            values << QString::number(1e10*Ascan_ds[k]);
                            if (std::isnan(1e10*Ascan_ds[k]))
                                qDebug() << "nan at: " << k;
                        }
                        out << values.join(',') << "\n";
                    }
                    this->progressBarPage3->setValue(100 * i / this->C_scan_double.size());
                    // Update the progress bar
                    QCoreApplication::processEvents(); // Allow GUI updates
                }
                // Reset or Update the progress bar for the next SNR value
                this->progressBarPage3->setValue(0);
                file.close();
            } else {
                // Handle invalid SNR input
                // You might want to add a message box or log here
            }
        }
        this->progressBarPage3->setValue(100); // Final progress bar update
    }
}

void ultrasound_cscan_seg::selectFolder() {
    QString dir = QFileDialog::getExistingDirectory(this, tr("Select Folder"), "/home",
                                                    QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    this->folderLabel->setText("Selected folder: " + dir);
    processFolder(dir);
}

void ultrasound_cscan_seg::processFolder(const QString &path) {
    QDir dir(path);
    // Get the list of all directories and files, excluding '.' and '..'
    QFileInfoList list = dir.entryInfoList(QDir::Files | QDir::Dirs | QDir::NoDotAndDotDot, QDir::DirsFirst);

    for (const QFileInfo &fileInfo : list) {
        if (fileInfo.isDir()) {
            // Process only the first-level subdirectories
            QDir subDir(fileInfo.absoluteFilePath());
            QFileInfoList subList = subDir.entryInfoList(QStringList() << "_*.csv", QDir::Files);
            for (const QFileInfo &subFileInfo : subList) {
                // Process .csv file and read the data
                qDebug() << "Found CSV file:" << subFileInfo.absoluteFilePath();
                this->fn = subFileInfo.absoluteFilePath();
                QFileInfo fileInfo(this->fn);
                this->C_scan_double.clear();
                this->C_scan_AS.clear();
                this->processFile(fileInfo);
                fillNanWithNearestNeighbors(this->C_scan_double);
                // align the surface
                ultrasound_Cscan_process::calculateSurface();
                // read downsample ratio
                bool ok_ds;
                int downsampleRate = this->downsampleRateInput->text().toInt(&ok_ds);
                if (!ok_ds || downsampleRate <= 0) {
                    // Handle error: invalid input
                    return;
                }
                int min_idx = 50*downsampleRate; // manual setting, times the ds factor
                ultrasound_Cscan_process::handleButton_alignsurface(min_idx);
                // save
                this->handleButton_multiSNR();
            }
        }
    }
}


// read the crop values (start, end) from the inputline
void ultrasound_cscan_seg::readCropValues(int& start, int& end) {
    // Default values
    start = 0;
    end = 0;
    if (this->cropSignalInput) {
        QString inputText = this->cropSignalInput->text();
        QStringList values = inputText.split(',');
        // Ensure two values are provided
        if (values.size() == 2) {
            bool ok1, ok2;
            // Parse the start value
            int tempStart = values[0].trimmed().toInt(&ok1);
            // Parse the end value
            int tempEnd = values[1].trimmed().toInt(&ok2);
            if (ok1 && ok2) {
                start = tempStart;
                end = tempEnd;
            } else {
                // Handle the error if the values are not integers
                qDebug() << "Invalid input for cropping values. Please enter integers.";
            }
        } else {
            // Handle the error if there aren't exactly two values
            qDebug() << "Please enter exactly two values for cropping, separated by a comma.";
        }
    } else {
        // Handle the error if cropSignalInput is not initialized
        qDebug() << "Crop Signal Input field is not initialized.";
    }
}
