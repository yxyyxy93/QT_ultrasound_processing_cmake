#include "ultrasound_cscan_process.h"
#include "qcustomplot.h"
#include "utils.h"
#include "trimdialog.h"
#include "OrthosliceViewer.h"

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
    , myButton_alignsurface(nullptr)
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
    connect(this->myButton_orthoslice, &QPushButton::clicked, this, [this]() {
        OrthosliceViewer *viewer = new OrthosliceViewer(nullptr,
                                                        this->C_scan_AS,
                                                        this->C_scan_double);
        viewer->setAttribute(Qt::WA_DeleteOnClose); // Ensure it's deleted on close
        viewer->show();
    });
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
        this->layout->addWidget(snrInput);
    }
    // Check if the Polynomial Degree Input Field already exists in the layout
    if (!widgetExistsInLayout<QLineEdit>(this->layout, "polynomialDegreeInput")) {
        // Create the Polynomial Degree Input Field
        QLineEdit* polynomialDegreeInput = new QLineEdit(this);
        polynomialDegreeInput->setObjectName("polynomialDegreeInput"); // Assign unique object name
        polynomialDegreeInput->setPlaceholderText("Enter degree of polynomial");
        this->layout->addWidget(polynomialDegreeInput);
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
        qDebug() << x << y << z;
        // Resize your member 3D QVector based on the dimensions
        this->C_scan_double.resize(x);
        for (int i = 0; i < x; ++i) {
            this->C_scan_double[i].resize(y);
            for (int j = 0; j < y; ++j) {
                QString line = in.readLine();
                QStringList values = line.split(',');
                // Resize the third dimension based on the actual number of values
                this->C_scan_double[i][j].resize(values.size());
                for (int index = 0; index < values.size(); ++index) {
                    this->C_scan_double[i][j][index] = values[index].toDouble();
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
    qDebug() << "finish trimming";
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
        // dynamic memory management
        this->pushButtons.append(myButton_plotsurface);
    }
    if (!this->myButton_alignsurface) { // Assuming myButton_alignsurface is declared in the header
        this->myButton_alignsurface = new QPushButton(tr("Align front surface"), this);
        this->layout->addWidget(myButton_alignsurface);
        int min_idx = 50; // manual setting
        connect(myButton_alignsurface,
                &QPushButton::clicked,
                this,
                [this, min_idx]() {handleButton_alignsurface(min_idx);}
                );
        this->pushButtons.append(myButton_alignsurface);
    }
}

void ultrasound_Cscan_process::handleButton_alignsurface(int min_idx){
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
    int threshold = 1000; // threshold for the max index
    for (int i = 0; i < this->C_scan_AS.size(); i++) {
        for (int j = 0; j < this->C_scan_AS[i].size(); j++) {
            QVector<std::complex<double>> Ascan_as = this->C_scan_AS[i][j];
            //            auto maxElementIndex = std::max_element(Ascan_as.begin(), Ascan_as.end(),
            //                                                    [](std::complex<double> a, std::complex<double> b) {
            //                                                        return std::abs(a) < std::abs(b);
            //                                                    });
            // *********** tailored for exp. *************** !
            double threshold_val = 0.5;
            auto maxElementIndex = std::find_if(Ascan_as.begin(), Ascan_as.end(),
                                                [threshold_val](std::complex<double> a) {
                                                    return std::abs(a) > threshold_val;
                                                });
            double currentIndex = std::distance(Ascan_as.begin(), maxElementIndex);
            if (currentIndex > threshold) {
                qDebug() << "out of threshold";
                currentIndex = prevMaxIndex;  // Use the previous value if threshold is exceeded
            } else {
                prevMaxIndex = currentIndex;  // Update the previous value
            }
            Front_surface_idx_i.push_back(currentIndex);
            Front_surface_val_i.push_back(std::abs(*maxElementIndex));
            // qDebug() << currentIndex;
        }
        this->Front_surface_idx.push_back(Front_surface_idx_i);
        this->Front_surface_val.push_back(Front_surface_val_i);
        Front_surface_idx_i.clear();
        Front_surface_val_i.clear();
    }
    // ******** 2D median filter or fitting
    QLineEdit* polynomialDegreeInput = this->findChild<QLineEdit*>("polynomialDegreeInput");
    int degree = 3;
    if (polynomialDegreeInput) {
        bool ok;
        int degree = polynomialDegreeInput->text().toInt(&ok);
        qDebug() << degree;
    }
    // take polynomialDegreeInput as filter size
    this->Front_surface_idx = applyMedianFilter(this->Front_surface_idx, degree);
    // this->Front_surface_idx = fitSurface(this->Front_surface_idx, degree);
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
        // for "if *==nullpter" checking
        this->myButton_alignsurface = nullptr;
        this->myButton_plotsurface = nullptr;
    }
    pushButtons.clear();

    for (QDoubleSpinBox* SpinBox : SpinBoxes) {
        delete SpinBox;
    }
    SpinBoxes.clear();
}

// // *********** for C-scan settings
// void ultrasound_Cscan_process::updateCscanPlotSelection(int index) {
//     this->CscanPlotMode = index; // This should be a member variable to store the current selection
//     qDebug() << index;
// }

