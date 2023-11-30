#include "ultrasound_cscan_process.h"
#include "qcustomplot.h"
#include "utils.h"
#include "trimdialog.h"

#include <QtWidgets>
#include <QHBoxLayout>
#include <QScrollBar>
#include <QDebug>
#include <QList>
#include <QTreeView>
#include <QStandardItemModel>
#include <QDir>
#include <QFile>
#include <algorithm>
#include <fstream>
#include <QFileDialog>

#include "npy.hpp"

ultrasound_Cscan_process::ultrasound_Cscan_process(QWidget *parent,
                                                   QString fn,
                                                   int fs,
                                                   double fx,
                                                   double fy)
    : QWidget(parent)
    , fn(fn)
    , fs(fs)
    , fx(fx)
    , fy(fy)
    , customPlot_Ascan(nullptr)
{
    this->layout = new QVBoxLayout(this);
    // Set the QVBoxLayout as the main layout of the widget
    this->setLayout(this->layout);

    // Ensure that the QWidget and its layout can expand vertically
    this->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->layout->setAlignment(Qt::AlignTop); // Align widgets to the top, allowing them to expand downwards

    // ************** create the push buttons and correponding labels
    this->myButton_load = new QPushButton(tr("Load data"), this);
    this->layout->addWidget(this->myButton_load);
    connect(this->myButton_load,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_load);
    //
    this->myButton_save = new QPushButton(tr("Save"), this);
    connect(this->myButton_save,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_save);
    this->layout->addWidget(this->myButton_save);
    //
    this->myButton_orthoslice = new QPushButton(tr("Orthoslice"), this);
    connect(this->myButton_orthoslice,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_orthoslice);
    this->layout->addWidget(this->myButton_orthoslice);
    //
    this->myButton_surface = new QPushButton("Determine the surface", this);
    connect(this->myButton_surface,
            &QPushButton::clicked,
            this,
            &ultrasound_Cscan_process::handleButton_surface);
    this->layout->addWidget(this->myButton_surface);

    // ************** progress bar
    this->m_progressBar = new QProgressBar;
    this->m_progressBar->setRange(0, 100);
    this->addNewWidgetAndReorderLayout(this->m_progressBar);

    // ****
}

ultrasound_Cscan_process::~ultrasound_Cscan_process(){
    // free any resources that were allocated by the class
    delete myButton_load;
    delete myButton_save;
    //    delete myButton_loadraw;
    delete myButton_orthoslice;
    delete myButton_surface;

    delete m_progressBar;

    delete myButton_alignsurface;
    delete myButton_plotsurface;

    delete customPlot1;
    delete customPlot2;
    delete customPlot3;
    delete scrollBarX;
    delete scrollBarY;
    delete scrollBarZ;

    delete customPlot_Ascan;

    delete customPlot_frontI;
    delete customPlot_frontV;

    delete layout;
}

// *********** read data - button
void ultrasound_Cscan_process::handleButton_load()
{
    // get the file name
    this->fn = QFileDialog::getOpenFileName(nullptr,
                                            "Open file",
                                            "",
                                            "Text files (*.txt; *.tdms; *.npy; *.bin; *.csv)");
    if (this->fn.isEmpty()) {
        qWarning() << "No file selected.";
    }
    // Determine the file type based on its extension
    QFileInfo fileInfo(this->fn);
    this->C_scan_double.clear();
    this->C_scan_AS.clear();
    this->processFile(fileInfo);
    qDebug() << this->C_scan_double.size();
    qDebug() << this->C_scan_double[0].size();
    qDebug() << this->C_scan_double[0][0].size();
    // fill nan
    fillNanWithNearestNeighbors(this->C_scan_double);
    QVector<std::complex<double>> Ascan_as;
    QVector<QVector<std::complex<double>>> Bscan_as;
    // calculate C_scan_AS
    for(int i = 0; i < this->C_scan_double.size(); i++) {
        for(int j = 0; j < this->C_scan_double[i].size(); j++) {
            QVector<double> Ascan = this->C_scan_double[i][j];
            Ascan_as = analyticSignal(Ascan);
            Bscan_as.push_back(Ascan_as);
        }
        this->m_progressBar->setValue(100 * i / this->C_scan_double.size());
        // Update the progress bar
        QCoreApplication::processEvents(); // Allow GUI updates
        this->C_scan_AS.push_back(Bscan_as);
        Bscan_as.clear();
    }
    this->m_progressBar->setValue(100);
    // add a button to trim the dataset
    if (!widgetExistsInLayout<QPushButton>(this->layout, "myButton_trim")) {
        QPushButton* myButton_trim = new QPushButton(tr("Trim the dataset"), this);
        myButton_trim->setObjectName("myButton_trim"); // assign unique object name
        this->layout->addWidget(myButton_trim);
        connect(myButton_trim,
                &QPushButton::clicked, this,
                &ultrasound_Cscan_process::handleButton_trim);
    }
    // Check if the 'Add Noise' button already exists in the layout
    if (!widgetExistsInLayout<QPushButton>(this->layout, "myButton_addNoise")) {
        // Create the 'Add Noise' button
        QPushButton* myButton_addNoise = new QPushButton(tr("Add Noise"), this);
        myButton_addNoise->setObjectName("myButton_addNoise"); // Assign unique object name
        this->layout->addWidget(myButton_addNoise);
        connect(myButton_addNoise, &QPushButton::clicked, this, &ultrasound_Cscan_process::handleButton_addNoise);
    }
    if (!widgetExistsInLayout<QLineEdit>(this->layout, "snrInput")) {
        QLineEdit* snrInput = new QLineEdit(this);
        snrInput->setObjectName("snrInput"); // Assign unique object name
        snrInput->setPlaceholderText("Enter SNR value");
        this->addNewWidgetAndReorderLayout(snrInput);
    }
}

void ultrasound_Cscan_process::handleButton_save(){
    QString filename = QFileDialog::getSaveFileName(this,
                                                    tr("Save Vector"),
                                                    QDir::homePath(),
                                                    tr("Text Files (*.txt)"));
    QFile file(filename);
    if (file.open(QIODevice::WriteOnly)) {
        QDataStream out(&file); // Create a QDataStream to write to the file

        // Write the dimensions as a header
        qint32 n1 = this->C_scan_double.size();
        qint32 n2 = this->C_scan_double[0].size();
        qint32 n3 = this->C_scan_double[0][0].size();
        out << n1 << n2 << n3;

        // Serialize the data to the stream
        out << this->C_scan_double;

        file.close();
    }
    // save C_scan_AS
    filename = QFileDialog::getSaveFileName(this,
                                            tr("Save Vector"),
                                            QDir::homePath(),
                                            tr("Text Files AS (*.txt)"));
    std::string stdString = filename.toStdString();
    std::ofstream fileas(stdString, std::ios::binary | std::ios::out);

    // Write the dimensions
    int32_t n1 = this->C_scan_AS.size();
    int32_t n2 = this->C_scan_AS[0].size();
    int32_t n3 = this->C_scan_AS[0][0].size();
    fileas.write(reinterpret_cast<char*>(&n1), sizeof(n1));
    fileas.write(reinterpret_cast<char*>(&n2), sizeof(n2));
    fileas.write(reinterpret_cast<char*>(&n3), sizeof(n3));

    // Write the data
    for (int i = 0; i < n1; ++i) {
        for (int j = 0; j < n2; ++j) {
            for (int k = 0; k < n3; ++k) {
                double value = abs(this->C_scan_AS[i][j][k]);
                fileas.write(reinterpret_cast<char*>(&value), sizeof(value));
            }
        }
    }
    file.close();
}

void ultrasound_Cscan_process::processFile(const QFileInfo &fileInfo) {
    QString extension = fileInfo.suffix();
    QVector<QVector<double>> B_scan_double;
    QVector<double> A_scan_double;
    if (extension == "txt") {
        // The file is a text file
        // read the data
        QFile file(this->fn);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qWarning() << "Could not open file:" << file.errorString();
        }
        QString contents = file.readAll().constData();
        file.close();
        // split the raw data
        QStringList contents_split;
        contents_split << contents.split("\n");
        // *************** read data
        QStringList lineData;
        // count the signals
        int count(1); // avoid divied by zero
        int lengthThreshold = 500;
        for (const QString &str : contents_split) {
            if (str.length() > lengthThreshold) {
                count++;
            }
        }
        qDebug() << "count of signals:" << count-1;
        //
        int count_cur(1);
        for (const QString &thisline: contents_split){
            // flag of the start of a B-scan
            if (thisline.length()>=5 && thisline.contains("GROUP")){
                C_scan_double.push_back(B_scan_double);
                B_scan_double.clear();
            }
            lineData.clear();
            lineData << thisline.split("\t");
            if (lineData.size() > lengthThreshold){
                count_cur++;
                m_progressBar->setValue(100 * count_cur / count);
                A_scan_double.clear();
                for (const QString &elem: lineData)
                {
                    A_scan_double.push_back(elem.toDouble());
                    // Update the progress bar
                    QCoreApplication::processEvents(); // Allow GUI updates
                }
                B_scan_double.push_back(A_scan_double);
            }
        }
        // last B_scan
        this->C_scan_double.removeFirst(); // remove the first one
        this->C_scan_double.push_back(B_scan_double);
        m_progressBar->setValue(100);
    } else if (extension == "csv"){
        QFile file(this->fn);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            qWarning() << "Could not open file:" << file.errorString();
        }
        QTextStream in(&file);
        // Read dimensions from the first line
        QString header = in.readLine();
        QStringList dimensions = header.split(',');
        int x = dimensions[0].toInt();
        int y = dimensions[1].toInt();
        int z = dimensions[2].toInt();
        // Resize your member 3D QVector based on the dimensions
        this->C_scan_double.resize(x);
        for (int i = 0; i < x; ++i) {
            this->C_scan_double[i].resize(y);
            for (int j = 0; j < y; ++j) {
                this->C_scan_double[i][j].resize(z);
            }
        }
        // Read the flattened data and reshape it
        for (int i = 0; i < x; ++i) {
            for (int j = 0; j < y; ++j) {
                QString line = in.readLine();
                QStringList values = line.split(',');
                for (int k = 0; k < z; ++k) {
                    this->C_scan_double[i][j][k] = values[k].toDouble();
                }
            }
        }
        file.close();
        m_progressBar->setValue(100);
    } else if (extension == "npy"){
        std::vector<unsigned long> shape {};
        bool fortran_order(0);
        std::vector<double> data;
        const std::string path = this->fn.toStdString();
        npy::LoadArrayFromNumpy(path, shape, fortran_order, data);
        //        double max_data = *max_element(data.begin(), data.end());
        // Define the dimensions of the output QVector !! to be improved !!!
        //        int x = 146; // for 2-ceh206-8-p21-0
        //        int y = 150;
        //        int rs_para = 0;
        int x = shape.at(0);
        int y = shape.at(1);
        int rs_para = 0;
        int z = data.size()/x/y;
        // Copy the input data into the output QVector
        int idx = 0;
        for (int i = 0; i < x; i++) {
            B_scan_double.clear();
            for (int j = 0; j < y; j++) {
                A_scan_double.clear();
                for (int k = 0; k < z; k++) {
                    idx++;
                    A_scan_double.push_back(data[idx]);
                }
                if(j>=rs_para && j<x+rs_para)
                    B_scan_double.push_back(A_scan_double);
            }
            m_progressBar->setValue(100 * i / x);
            // Update the progress bar
            QCoreApplication::processEvents(); // Allow GUI updates
            this->C_scan_double.push_back(B_scan_double);
        }
        //        // compensate to make it square
        //        for (int i = 0; i< y-x; i++)
        //            this->C_scan_double.push_back(B_scan_double);
        m_progressBar->setValue(100);
    } else if(extension == "bin"){ // "Load .bin files from NTU"
        QString path = fileInfo.absolutePath();
        QDir dir(path);
        if (!dir.exists()) {
            qWarning() << "Directory does not exist:" << path;
            return;
        }
        // read the .dat file
        dir.setNameFilters(QStringList() << "*.dat");
        dir.setFilter(QDir::Files | QDir::NoDotAndDotDot | QDir::NoSymLinks);
        QFileInfoList fileInfoList = dir.entryInfoList();
        int tl(0);
        if(fileInfoList.size() == 1) {
            QVector<double> data = readDatFile(fileInfoList.first().absoluteFilePath());
            // Get the length of the data
            tl = data.size();
            // Display the length
            qDebug() << "Length of data:" << tl;
        } else if(fileInfoList.size() > 1) {
            qWarning() << "More than one .dat file in directory.";
        } else {
            qWarning() << "No .dat file found in directory.";
        }
        // read the .bin files
        dir.setNameFilters(QStringList() << "*.bin");
        dir.setFilter(QDir::Files | QDir::NoDotAndDotDot | QDir::NoSymLinks);
        fileInfoList = dir.entryInfoList();
        // Sorting files by creation time
        std::sort(fileInfoList.begin(), fileInfoList.end(), [](const QFileInfo &file1, const QFileInfo &file2){
            return file1.lastModified() < file2.lastModified();
        });
        for (const QFileInfo& fileInfo : fileInfoList) {
            QFile file(fileInfo.absoluteFilePath());
            if (!file.open(QIODevice::ReadOnly)) {
                qWarning() << "Cannot open file:" << fileInfo.fileName() << " Error:" << file.errorString();
                continue;
            }
            // Reading data from file
            QDataStream in(&file);
            // Set the byte order to LittleEndian if the file was written on a little-endian machine
            in.setByteOrder(QDataStream::LittleEndian);
            Array1D data;
            while (!in.atEnd()) {
                double doubleValue;
                in >> doubleValue;
                data.append(doubleValue);
            }
            file.close();
            // Processing data
            // reshape the data
            Array2D reshaped(data.size()/tl, Array1D(tl, 0));
            for (int i = 0; i < data.size()/tl; ++i) {
                for (int j = 0; j < tl; ++j) {
                    reshaped[i][j] = data[i*tl + j];
                }
            }
            this->C_scan_double.push_back(reshaped);
            qDebug() << "Read file:" << fileInfo.fileName() << "Size:" << data.size();
        }
    } else {
        qDebug() << "The file type is unknown";
    }
}

// manipulate dataset
void ultrasound_Cscan_process::handleButton_trim(){
    TrimDialog dialog(this);
    if (dialog.exec() == QDialog::Accepted) {
        int startI = dialog.startI();
        int endI = dialog.endI();
        int startJ = dialog.startJ();
        int endJ = dialog.endJ();
        int startK = dialog.startK();
        int endK = dialog.endK();

        // Now you can use the input values to trim your 3D QVector...
        this->C_scan_AS = trim3DData(this->C_scan_AS,
                                     startI, endI,
                                     startJ, endJ,
                                     startK, endK);
        this->C_scan_double = trim3DData(this->C_scan_double,
                                         startI, endI,
                                         startJ, endJ,
                                         startK, endK);
    }
}

void ultrasound_Cscan_process::handleButton_addNoise(){
    QLineEdit* snrInput = this->findChild<QLineEdit*>("snrInput");
    if (snrInput) {
        bool ok;
        double snrDb = snrInput->text().toDouble(&ok);
        if (ok) {
            QVector<std::complex<double>> Ascan_as;
            // calculate C_scan_AS
            for(int i = 0; i < this->C_scan_double.size(); i++) {
                for(int j = 0; j < this->C_scan_double[i].size(); j++) {
                    QVector<double> Ascan = this->C_scan_double[i][j];
                    addGaussianNoise(Ascan, snrDb);
                    Ascan_as = analyticSignal(Ascan);
                    this->C_scan_AS[i][j] = Ascan_as;
                }
                this->m_progressBar->setValue(100 * i / this->C_scan_double.size());
                // Update the progress bar
                QCoreApplication::processEvents(); // Allow GUI updates
            }
        } else {
            // Handle invalid SNR input
        }
        m_progressBar->setValue(100);
    }
}

// **************** define surfac
void ultrasound_Cscan_process::handleButton_surface()
{
    this->calculateSurface();

    this->m_progressBar->setValue(100);
    // add a Qlabel
    this->myButton_surface ->setText("Determine the surface (Surface found!)");
    // add a button to plot the surface
    // ************** create the push buttons and correponding labels
    // Check if the buttons already exist before creating them
    if (!this->myButton_plotsurface) { // Assuming myButton_plotsurface is declared in the header
        this->myButton_plotsurface = new QPushButton(tr("Plot front surface"), this);
        this->layout->addWidget(myButton_plotsurface);
        connect(this->myButton_plotsurface,
                &QPushButton::clicked, this,
                &ultrasound_Cscan_process::handleButton_plotsurface);
    }
    if (!this->myButton_alignsurface) { // Assuming myButton_alignsurface is declared in the header
        this->myButton_alignsurface = new QPushButton(tr("Align front surface"), this);
        this->layout->addWidget(myButton_alignsurface);
        connect(myButton_alignsurface,
                &QPushButton::clicked, this,
                &ultrasound_Cscan_process::handleButton_alignsurface);
    }

    // dynamic memory management
    this->pushButtons.append(myButton_plotsurface);
    this->pushButtons.append(myButton_alignsurface);
}

void ultrasound_Cscan_process::handleButton_alignsurface(){
    // get min of the front surface index
    int min_idx = this->C_scan_AS[0][0].size();
    for(int i = 0; i < this->Front_surface_idx.size(); i++) {
        for(int j = 0; j < this->Front_surface_idx[i].size(); j++) {
            min_idx = (min_idx>=this->Front_surface_idx[i][j])?
                          this->Front_surface_idx[i][j]:min_idx;
        }
    }
    min_idx = 500; // manual setting !!!!!!!!
    qDebug() << min_idx;
    // shift
    for(int i = 0; i < this->Front_surface_idx.size(); i++){
        for(int j = 0; j < this->Front_surface_idx[i].size(); j++){
            int front_idx = this->Front_surface_idx[i][j];
            // qDebug() << front_idx;
            shiftVector_1D(this->C_scan_AS[i][j],
                           front_idx-min_idx);
            shiftVector_1D(this->C_scan_double[i][j],
                           front_idx-min_idx);
        }
    }
    if(this->myButton_alignsurface)
        this->myButton_alignsurface ->setText("Align front surface (Surface aligned!)");
}

void ultrasound_Cscan_process::handleButton_plotsurface(){
    // plot the surface
    QCustomPlot *customPlot_fsurface_idx = new QCustomPlot();
    QCustomPlot *customPlot_fsurface_val = new QCustomPlot();
    //
    this->layout->addWidget(customPlot_fsurface_idx);
    this->layout->addWidget(customPlot_fsurface_val);
    // Create QCPColorMap objects
    QCPColorMap *map1 = new QCPColorMap(customPlot_fsurface_idx->xAxis,
                                        customPlot_fsurface_idx->yAxis);
    QCPColorMap *map2 = new QCPColorMap(customPlot_fsurface_val->xAxis,
                                        customPlot_fsurface_val->yAxis);
    // set the data for each QCPColorMap
    int x_size = this->Front_surface_idx.size();
    int y_size = this->Front_surface_idx[0].size();
    //
    map1->data()->setSize(x_size, y_size);
    map1->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    map2->data()->setSize(x_size, y_size);
    map2->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    //
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            map1->data()->setCell(i, j, this->Front_surface_idx[i][j]);
        }
    }
    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            map2->data()->setCell(i, j, this->Front_surface_val[i][j]);
        }
    }
    map1->setGradient(QCPColorGradient::gpHot);
    map2->setGradient(QCPColorGradient::gpHot);
    map1->setInterpolate(true);
    map2->setInterpolate(true);
    // Add a color scale to the custom plot widget
    QCPColorScale *colorScale;
    colorScale = new QCPColorScale(customPlot_fsurface_idx);
    map1->setColorScale(colorScale);
    colorScale = new QCPColorScale(customPlot_fsurface_val);
    map2->setColorScale(colorScale);
    // Rescale the color map data range to fit the new data
    map1->rescaleDataRange();
    map2->rescaleDataRange();
    //
    customPlot_fsurface_idx->xAxis->setRange(0, x_size);
    customPlot_fsurface_idx->yAxis->setRange(0, y_size);
    customPlot_fsurface_val->xAxis->setRange(0, x_size);
    customPlot_fsurface_val->yAxis->setRange(0, y_size);
    // add adjustable colorscale
    QDoubleSpinBox *minColorBound = new QDoubleSpinBox();
    QDoubleSpinBox *maxColorBound = new QDoubleSpinBox();
    minColorBound->setRange(-1e9, 1e9); // Set appropriate limits
    maxColorBound->setRange(-1e9, 1e9); // Set appropriate limits
    this->layout->addWidget(minColorBound);
    this->layout->addWidget(maxColorBound);
    QObject::connect(minColorBound, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                     [map1](double newMin) {
                         map1->setDataRange(QCPRange(newMin, map1->dataRange().upper));
                         map1->parentPlot()->replot();
                     });
    QObject::connect(maxColorBound, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                     [map1](double newMax) {
                         map1->setDataRange(QCPRange(map1->dataRange().lower, newMax));
                         map1->parentPlot()->replot();
                     });
    // Call replot() to update the plot with the new data
    customPlot_fsurface_idx->replot();
    customPlot_fsurface_val->replot();

    // dynamic memory management
    this->customPlots.append(customPlot_fsurface_idx);
    this->customPlots.append(customPlot_fsurface_val);

    this->maps.append(map1);
    this->maps.append(map2);

    this->colorScales.append(colorScale);

    this->SpinBoxes.append(minColorBound);
    this->SpinBoxes.append(maxColorBound);
    // ...
}

void ultrasound_Cscan_process::calculateSurface() {
    if (!this->Front_surface_idx.isEmpty())
        this->Front_surface_idx.clear();
    if (!this->Front_surface_val.isEmpty())
        this->Front_surface_val.clear();

    QVector<double> Front_surface_idx_i;
    QVector<double> Front_surface_val_i;
    // I set this threshold to avoid the bug in simulation data: finding the second echo as max.
    int prevMaxIndex = 0;  // Initial value, adjust based on your data requirements
    int threshold = 600; // threshold for the max index
    for (int i = 0; i < this->C_scan_AS.size(); i++) {
        for (int j = 0; j < this->C_scan_AS[i].size(); j++) {
            QVector<std::complex<double>> Ascan_as = this->C_scan_AS[i][j];
            auto maxElementIndex = std::max_element(Ascan_as.begin(), Ascan_as.end(),
                                                    [](std::complex<double> a, std::complex<double> b) {
                                                        return std::abs(a) < std::abs(b);
                                                    });

            int currentIndex = std::distance(Ascan_as.begin(), maxElementIndex);
            if (currentIndex > threshold) {
                currentIndex = prevMaxIndex;  // Use the previous value if threshold is exceeded
            } else {
                prevMaxIndex = currentIndex;  // Update the previous value
            }
            Front_surface_idx_i.push_back(currentIndex);
            Front_surface_val_i.push_back(std::abs(*maxElementIndex));
        }
        this->Front_surface_idx.push_back(Front_surface_idx_i);
        this->Front_surface_val.push_back(Front_surface_val_i);
        Front_surface_idx_i.clear();
        Front_surface_val_i.clear();
    }
}

// *************visualization
void ultrasound_Cscan_process::handleButton_orthoslice() {
    // ****************** create the orthoslice visual
    // Create the QCustomPlot widget
    this->customPlot1 = new QCustomPlot();
    this->customPlot2 = new QCustomPlot();
    this->customPlot3 = new QCustomPlot();

    // Create the ScrollBars and Labels
    this->scrollBarX = new QScrollBar(Qt::Horizontal);
    this->scrollBarY = new QScrollBar(Qt::Horizontal);
    this->scrollBarZ = new QScrollBar(Qt::Horizontal);

    QLabel *sBX_label = new QLabel();
    QLabel *sBY_label = new QLabel();
    QLabel *sBZ_label = new QLabel();

    // ... [Connect signals and slots as before] ...
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarX, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarX->minimum()).arg(scrollBarX->maximum()).arg(value);
        sBX_label->setText(labelText);
    });
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarY, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarY->minimum()).arg(scrollBarY->maximum()).arg(value);
        sBY_label->setText(labelText);
    });
    // Connect the scrollbar's valueChanged() signal to a slot that updates the label text
    connect(scrollBarZ, &QScrollBar::valueChanged, [=](int value) {
        QString labelText = QString("Min: %1 Max: %2 Current: %3").
                            arg(scrollBarZ->minimum()).arg(scrollBarZ->maximum()).arg(value);
        sBZ_label->setText(labelText);
    });
    // Create the combo box and add the selection options
    QComboBox *comboBox = new QComboBox();
    comboBox->addItem("Max");
    comboBox->addItem("Average");
    comboBox->addItem("Single Index");
    comboBox->addItem("Single Index for original");
    // Connect the currentIndexChanged signal to a slot
    connect(comboBox, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, &ultrasound_Cscan_process::updateCscanPlotSelection);
    //
    QPushButton *deleteAllButton = new QPushButton("Delete All");
    this->addNewWidgetAndReorderLayout(deleteAllButton);
    connect(deleteAllButton, &QPushButton::clicked, this, &ultrasound_Cscan_process::clearAllDynamicMemory);
    // Connect the valueChanged() signals of the QScrollBars to update the plot data
    QObject::connect(this->scrollBarX, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);
    QObject::connect(this->scrollBarY, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);
    QObject::connect(this->scrollBarZ, &QScrollBar::valueChanged, this,
                     &ultrasound_Cscan_process::updatePlot);

    // Create the main horizontal layout
    QHBoxLayout *hLayout = new QHBoxLayout();
    // Assuming 'this->layout' is the existing layout of the widget
    this->layout->addLayout(hLayout);

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
    plot3l->addWidget(this->scrollBarZ);
    plot3l->addWidget(sBZ_label);
    plot3l->addWidget(comboBox); // Assuming comboBox is defined earlier
    plot3l->addWidget(this->customPlot3);

    // Add the vertical layouts to the main horizontal layout
    hLayout->addLayout(plot1l);
    hLayout->addLayout(plot2l);
    hLayout->addLayout(plot3l);
    // Add the horizontal layout to the existing main layout of the widget
    this->setLayout(hLayout);

    // Set the range of the QScrollBars based on the size of the data
    this->scrollBarX->setRange(0, this->C_scan_AS.size() - 1);
    this->scrollBarY->setRange(0, this->C_scan_AS[0].size() - 1);
    this->scrollBarZ->setRange(0, this->C_scan_AS[0][0].size() - 1);
    // Set the initial values of the QScrollBars
    this->scrollBarX->setValue(this->C_scan_AS.size() / 2);
    this->scrollBarY->setValue(this->C_scan_AS[0].size() / 2);
    this->scrollBarZ->setValue(this->C_scan_AS[0][0].size() / 2);
    // Set size policies for custom plots
    this->customPlot1->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->customPlot2->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->customPlot3->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    this->scrollBarX->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    this->scrollBarY->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    this->scrollBarZ->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sBX_label->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sBY_label->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    sBZ_label->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
    comboBox->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);

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

    plot2l->addWidget(scrollBarY);
    plot2l->addWidget(sBY_label);
    plot2l->addWidget(this->customPlot2);

    plot3l->addWidget(scrollBarZ);
    plot3l->addWidget(sBZ_label);
    plot3l->addWidget(comboBox);
    plot3l->addWidget(this->customPlot3);

    // Adjust the horizontal layout
    hLayout->setStretch(0, 1); // Assign stretch factor to plot1w
    hLayout->setStretch(1, 1); // Assign stretch factor to plot2w
    hLayout->setStretch(2, 1); // Assign stretch factor to plot3w

    // Set the horizontal layout to the main widget or window
    //    printWidgetInfo(this->customPlot3);
    //    printLayoutInfo(plot3l);

    qDebug() << "Plot Layout size hint:" << plot1l->sizeHint();
    qDebug() << "CustomPlot1 size policy:" << this->customPlot1->sizePolicy().horizontalPolicy()
             << this->customPlot1->sizePolicy().verticalPolicy();

    // After making changes to the layout, update the widget
    this->update(); // This will schedule a repaint event for the widget
    this->adjustSize(); // This will adjust the size of the widget to fit its contents

    // dynamic memory management
    this->customPlots.append(this->customPlot1);
    this->customPlots.append(this->customPlot2);
    this->customPlots.append(this->customPlot3);

    this->scrollBars.append(this->scrollBarX);
    this->scrollBars.append(this->scrollBarY);
    this->scrollBars.append(this->scrollBarZ);

    this->labels.append(sBX_label);
    this->labels.append(sBY_label);
    this->labels.append(sBZ_label);

    this->hLayouts.append(hLayout);
    this->vLayouts.append(plot1l);
    this->vLayouts.append(plot2l);
    this->vLayouts.append(plot3l);

    this->pushButtons.append(deleteAllButton);

    this->comboxes.append(comboBox);
    // ...
}

void ultrasound_Cscan_process::updatePlot() {
    int sliceX = this->scrollBarX->value();
    int sliceY = this->scrollBarY->value();
    int sliceZ = this->scrollBarZ->value();
    qDebug() << sliceX;
    qDebug() << sliceY;
    qDebug() << sliceZ;
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
    //
    map1->data()->setSize(y_size, z_size);
    map1->data()->setRange(QCPRange(0, y_size),
                           QCPRange(0, z_size));
    map2->data()->setSize(x_size, z_size);
    map2->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, z_size));
    map3->data()->setSize(x_size, y_size);
    map3->data()->setRange(QCPRange(0, x_size),
                           QCPRange(0, y_size));
    //
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
                for (int idx = std::max(sliceZ - 20, 0); idx < std::min(sliceZ + 20, z_size); ++idx) {
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
    //
    map3->setGradient(QCPColorGradient::gpJet);
//    // Add a color scale to the custom plot widget
//    QCPColorScale *colorScale = new QCPColorScale(this->customPlot3);
//    // Set the color map for the color scale
//    colorScale->setDataRange(map3->dataRange());
//    colorScale->setGradient(map3->gradient());
//    //    map1->setColorScale(colorScale);
//    //    map2->setColorScale(colorScale);
//    map3->setColorScale(colorScale);
//    // add a color scale:
//    // Check if a color scale already exists
//    for (int i = 0; i < this->customPlot3->plotLayout()->elementCount(); ++i) {
//        QCPLayoutElement *element = this->customPlot3->plotLayout()->elementAt(i);
//        if (QCPColorScale *existingColorScale = qobject_cast<QCPColorScale *>(element)) {
//            // Remove existing color scale
//            this->customPlot3->plotLayout()->remove(existingColorScale);
//            break; // Assuming there is only one color scale
//        }
//    }
//    this->customPlot3->plotLayout()->addElement(0, 1, colorScale); // add it to the right of the main axis rect
//    colorScale->setType(QCPAxis::atRight); // scale shall be vertical bar with tick/axis labels right (actually atRight is already the default)// associate the color map with the color scale
//    colorScale->axis()->setLabel("Amp. (arb.)");

    // Rescale the color map data range to fit the new data
    map1->rescaleDataRange(true);
    map2->rescaleDataRange(true);
    map3->rescaleDataRange(true);
    //
    map1->rescaleAxes();
    map2->rescaleAxes();
    map3->rescaleAxes();
    //
    this->customPlot1->xAxis->setRange(0, y_size);
    this->customPlot1->yAxis->setRange(0, z_size);
    this->customPlot2->xAxis->setRange(0, x_size);
    this->customPlot2->yAxis->setRange(0, z_size);
    this->customPlot3->xAxis->setRange(0, x_size);
    this->customPlot3->yAxis->setRange(0, y_size);
    //    // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
    //    this->customPlot3->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    // Call replot() to update the plot with the new data
    this->customPlot1->replot();
    this->customPlot2->replot();
    this->customPlot3->replot();
    // Connect the plot's mouse press signal to the onCustomPlotClicked slot
    connect(this->customPlot3, &QCustomPlot::mousePress,
            this, &ultrasound_Cscan_process::onCustomPlotClicked_Cscan);
}

void ultrasound_Cscan_process::onCustomPlotClicked_Cscan(QMouseEvent* event)
{
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
        QVector<double> signal = this->C_scan_double[x][y];
        QVector<double> time;
        for (int i = 0; i < signal.size(); ++i) {
            time.append(i);
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
        if (this->customPlot_Ascan==nullptr){
            this->customPlot_Ascan = new QCustomPlot();
            this->layout->addWidget(this->customPlot_Ascan);
        }
        else
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
        this->customPlot_Ascan->graph(0)->setData(time,
                                                  signal);
        this->customPlot_Ascan->graph(1)->setData(time,
                                                  Ascan_as_abs);
        this->customPlot_Ascan->rescaleAxes();

        // Allow user to drag axis ranges with mouse, zoom with mouse wheel and select graphs by clicking:
        this->customPlot_Ascan->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
        this->customPlot_Ascan->replot();
    }
    // dynamic memory management
    this->customPlots.append(this->customPlot_Ascan);
}

void ultrasound_Cscan_process::addNewWidgetAndReorderLayout(QWidget* newWidget) {
    // Check if the newWidget already exists
    for (int i = 0; i < this->layout->count(); ++i) {
        QWidget* existingWidget = this->layout->itemAt(i)->widget();
        if (existingWidget == newWidget) {
            // Widget already exists in the layout, so don't add it again
            return;
        }
    }
    // Add the new widget
    this->layout->addWidget(newWidget);
    // Check if the vertical spacer already exists
    bool spacerFound = false;
    for (int i = 0; i < this->layout->count(); ++i) {
        if (this->layout->itemAt(i)->spacerItem()) {
            spacerFound = true;
            // Move the spacer to the end of the layout
            this->layout->addItem(this->layout->takeAt(i));
            break;
        }
    }
    // If no spacer found, add a new one
    if (!spacerFound) {
        QSpacerItem* verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Expanding, QSizePolicy::Expanding);
        this->layout->addItem(verticalSpacer);
    }
}

// ************** clear dynamic memory allocations
void ultrasound_Cscan_process::clearAllDynamicMemory() {
    // you have to delete colorscale, map, customplot in sequence. otherwise it crushes
    for (QCPColorScale *colorScale : colorScales) {
        delete colorScale;
    }
    colorScales.clear();

    for (QCPColorMap* map : maps) {
        map->deleteLater();
    }
    maps.clear();

    for (QCustomPlot *plot : customPlots) {
        delete plot;
    }
    customPlots.clear();

    for (QScrollBar *scrollBar : scrollBars) {
        delete scrollBar;
    }
    scrollBars.clear();

    for (QLabel *label : labels) {
        delete label;
    }
    labels.clear();

    for (QComboBox * combox: comboxes){
        delete combox;
    }
    comboxes.clear();

    for (QWidget *widget : widgets) {
        widget->deleteLater();
    }
    widgets.clear();

    for (QPushButton *pushButton : pushButtons) {
        delete pushButton;
    }
    pushButtons.clear();

    for (QDoubleSpinBox* SpinBox : SpinBoxes) {
        delete SpinBox;
    }
    SpinBoxes.clear();
}

// *********** for C-scan settings
void ultrasound_Cscan_process::updateCscanPlotSelection(int index) {
    this->CscanPlotMode = index; // This should be a member variable to store the current selection
    qDebug() << index;
}

