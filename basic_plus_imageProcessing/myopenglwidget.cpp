#include "MyOpenGLWidget.h"
#include <QMouseEvent>
#include <QScrollBar>
#include <QVBoxLayout>
#include <QPainter>
#include <QCheckBox>

#include "..\basic_read\utils.h"
#include "graphwindow.h"

// Adjust these values to position the axes at a desired location in the scene
float axisOffsetX = 1.5f;
float axisOffsetY = 1.5f;
float axisOffsetZ = 1.5f;

// Example axis vertices with offset
const GLfloat axisVertices[] = {
    // X axis (Red)
    axisOffsetX - 1.0f, axisOffsetY, axisOffsetZ, 1.0f, 0.0f, 0.0f,
    axisOffsetX + 1.0f, axisOffsetY, axisOffsetZ, 1.0f, 0.0f, 0.0f,
    // Y axis (Green)
    axisOffsetX, axisOffsetY - 1.0f, axisOffsetZ, 0.0f, 1.0f, 0.0f,
    axisOffsetX, axisOffsetY + 1.0f, axisOffsetZ, 0.0f, 1.0f, 0.0f,
    // Z axis (Blue)
    axisOffsetX, axisOffsetY, axisOffsetZ - 1.0f, 0.0f, 0.0f, 1.0f,
    axisOffsetX, axisOffsetY, axisOffsetZ + 1.0f, 0.0f, 0.0f, 1.0f
};
GLuint axisVAO, axisVBO;

// Constructor
MyOpenGLWidget::MyOpenGLWidget(QWidget *parent)
    : QOpenGLWidget(parent), m_program(nullptr), viewAngle_x(45.0f), viewAngle_y(45.0f), visibilityThreshold(0.5f), useDenoisedData(false) {

    // Initialize the scrollbar for adjusting the view angle
    QScrollBar *angleScrollBar_x = new QScrollBar(Qt::Horizontal, this);
    angleScrollBar_x->setRange(0, 360); // Range of angles
    angleScrollBar_x->setValue(45); // Initial angle
    QScrollBar *angleScrollBar_y = new QScrollBar(Qt::Horizontal, this);
    angleScrollBar_y->setRange(0, 360); // Range of angles
    angleScrollBar_y->setValue(45); // Initial angle

    // Connect the scrollbar's signal to the slot for angle change
    connect(angleScrollBar_x, &QScrollBar::valueChanged, this, &MyOpenGLWidget::setViewAngle_x);
    connect(angleScrollBar_y, &QScrollBar::valueChanged, this, &MyOpenGLWidget::setViewAngle_y);

    // Layout to add scrollbar below the OpenGL widget
    QVBoxLayout *layout = new QVBoxLayout(this);
    layout->addWidget(angleScrollBar_x);
    layout->addWidget(angleScrollBar_y);
    layout->setAlignment(Qt::AlignBottom);

    // Initialize the scrollbar for setting the visibility threshold
    QScrollBar *thresholdScrollBar = new QScrollBar(Qt::Horizontal, this);
    thresholdScrollBar->setRange(0, 100); // Threshold range (e.g., 0-100%)
    thresholdScrollBar->setValue(100); // Initial threshold value
    // Connect the scrollbar's signal to the slot for threshold change
    connect(thresholdScrollBar, &QScrollBar::valueChanged, this, &MyOpenGLWidget::setVisibilityThreshold);

    // Create and add the checkbox for useDenoisedData
    QCheckBox *denoiseCheckbox = new QCheckBox("Use Denoised Data", this);
    connect(denoiseCheckbox, &QCheckBox::stateChanged, this, &MyOpenGLWidget::toggleDenoisedData);
    layout->addWidget(denoiseCheckbox);

    // Set the layout alignment
    layout->setAlignment(Qt::AlignBottom);

    layout->addWidget(thresholdScrollBar);
    layout->setAlignment(Qt::AlignBottom);
}

// Slot to toggle useDenoisedData
void MyOpenGLWidget::toggleDenoisedData(int state) {
    useDenoisedData = (state == Qt::Checked);
    update(); // Trigger a repaint to reflect the change
}

// Slot to update the view angle
void MyOpenGLWidget::setViewAngle_x(int angle) {
    this->viewAngle_x = angle;
    update(); // Trigger a repaint
}
void MyOpenGLWidget::setViewAngle_y(int angle) {
    this->viewAngle_y = angle;
    update(); // Trigger a repaint
}

// Slot to update the visibility threshold
void MyOpenGLWidget::setVisibilityThreshold(int value) {
    this->visibilityThreshold = static_cast<float>(value) / 100.0f; // Convert to a percentage
    update(); // Trigger a repaint
}

MyOpenGLWidget::~MyOpenGLWidget() {
    makeCurrent();
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    delete m_program;
    doneCurrent();
}

void MyOpenGLWidget::initializeGL() {
    initializeOpenGLFunctions();
    // Set up the rendering context, define display lists etc.:
    static const char *vertexShaderSource =
        "#version 330 core\n"
        "layout (location = 0) in vec3 aPos;\n"
        "layout (location = 1) in vec3 aColor;\n" // input layout for vertex color
        "uniform mat4 model;\n"
        "uniform mat4 view;\n"
        "uniform mat4 projection;\n"
        "out vec3 ourColor;\n" // output for color to the fragment shader
        "void main() {\n"
        "   gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
        "   ourColor = aColor;\n" // Pass color to the fragment shader
        "}\0";

    static const char *fragmentShaderSource =
        "#version 330 core\n"
        "in vec3 ourColor;\n" // input from vertex shader
        "out vec4 FragColor;\n"
        "void main() {\n"
        "   FragColor = vec4(ourColor, 1.0);\n" //  Fragment color from vertex shader
        "}\n\0";

    m_program = new QOpenGLShaderProgram(this);
    m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
    m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);

    if (!m_program->link()) {
        qDebug() << "Shader Program Error:" << m_program->log();
        return; // Add error handling
    }

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, this->vertexData.size() * sizeof(GLfloat), this->vertexData.constData(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindBuffer(GL_ARRAY_BUFFER, 0); // Unbind the VBO

    // *************** for axis coordinates
    // Generate and bind the VAO for the axes
    glGenVertexArrays(1, &axisVAO);
    glBindVertexArray(axisVAO);
    // Generate and bind the VBO for the axes
    glGenBuffers(1, &axisVBO);
    glBindBuffer(GL_ARRAY_BUFFER, axisVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(axisVertices), axisVertices, GL_STATIC_DRAW);
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)0);
    glEnableVertexAttribArray(0);
    // Color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (GLvoid*)(3 * sizeof(GLfloat)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0); // Unbind VAO

    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
}

void MyOpenGLWidget::resizeGL(int w, int h) {
    // Calculate aspect ratio
    qreal aspect = qreal(w) / qreal(h ? h : 1);
    // Set near plane to 0.1, far plane to 100.0, field of view 45 degrees
    const qreal zNear = 0.1, zFar = 100.0, fov = 45.0;

    // Reset projection
    QMatrix4x4 projection;
    projection.perspective(fov, aspect, zNear, zFar);
    m_program->bind(); // Bind to set uniforms
    m_program->setUniformValue("projection", projection);
    m_program->release(); // Release after setting uniforms
}

void MyOpenGLWidget::paintGL() {
    // Clear the buffer
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set clear color to white
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
    glPointSize(5.0f); // Set the point size

    m_program->bind(); // Bind your shader program
    // Set up transformation matrices
    QMatrix4x4 model;
    model.translate(0.0f, 0.0f, -2.0f);
    model.rotate(this->viewAngle_x, 0.0f, 1.0f, 0.0f); // Rotate around the y-axis using the scrollbar angle
    model.rotate(this->viewAngle_y, 1.0f, 0.0f, 0.0f); // Rotate around the y-axis using the scrollbar angle

    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f);

    m_program->setUniformValue("model", model);
    m_program->setUniformValue("view", view);
    // Choose the data set based on the useDenoisedData flag
    const QVector<GLfloat>& currentVertexData = this->useDenoisedData ? vertexData_denoise : vertexData;

    // Draw the point cloud or other objects
    glBindVertexArray(VAO);
    int numPoints = currentVertexData.size() / 6; // Assuming 6 floats per vertex
    for (int i = 0; i < numPoints; i++) {
        int vertexStart = i * 6;
        if (currentVertexData[vertexStart + 3] >= this->visibilityThreshold) {
            glDrawArrays(GL_POINTS, i, 1); // Draw each point that meets the threshold
        }
    }
    glBindVertexArray(0); // Unbind the VAO
    // Draw the axes
    glBindVertexArray(axisVAO); // Ensure you have created this VAO and uploaded the axis vertices
    glDrawArrays(GL_LINES, 0, 6); // 6 vertices total for the 3 lines (2 per line)
    glBindVertexArray(0); // Unbind the axis VAO
    m_program->release(); // Release the shader program
}

// ************* convert the structure to vertices
void MyOpenGLWidget::convertStructureToVertices(const QVector<QVector<QVector<double>>>& structure) {
    QVector<GLfloat> vertices;
    // Get the maximum size for normalization
    int sizeX = structure.size();
    int sizeY = sizeX > 0 ? structure[0].size() : 0;
    int sizeZ = (sizeX > 0 && sizeY > 0) ? structure[0][0].size() : 0;
    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();
    // Find the min and max values in the dataset
    for (const auto& layer : structure) {
        for (const auto& row : layer) {
            for (const auto& val : row) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
            }
        }
    }
    int step = 3; // Step size of 3
    // Now, minVal and maxVal can be used in the mapValueToColor function
    for (int x = 0; x < sizeX; x += step) { // Increment x by step
        for (int y = 0; y < sizeY; y += step) { // Increment y by step
            for (int z = 0; z < sizeZ; z += step) { // Increment z by step
                // Normalize the indices to the range [-1, 1]
                // Adjust calculation to account for downsampling
                float posX = 2.0f * x / (sizeX - 1) - 1.0f;
                float posY = 2.0f * y / (sizeY - 1) - 1.0f;
                float posZ = 2.0f * z / (sizeZ - 1) - 1.0f;
                // Add position data
                this->vertexData.push_back(posX);
                this->vertexData.push_back(posY);
                this->vertexData.push_back(posZ);
                // Assuming mapValueToColor is a function that maps a data value to a color
                // And assuming structure[x][y][z] is valid and contains the data value you want to map
                QVector3D color = mapValueToColor(structure[x][y][z], minVal, maxVal);
                this->vertexData.push_back(color.x()); // R
                this->vertexData.push_back(color.y()); // G
                this->vertexData.push_back(color.z()); // B
            }
        }
    }
    QVector<QVector<QVector<double>>> structure_dn = denoise3D(structure, 1);
    // Find the min and max values in the dataset
    for (const auto& layer : structure_dn) {
        for (const auto& row : layer) {
            for (const auto& val : row) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
            }
        }
    }
    // Now, minVal and maxVal can be used in the mapValueToColor function
    for (int x = 0; x < sizeX; x += step) { // Increment x by step
        for (int y = 0; y < sizeY; y += step) { // Increment y by step
            for (int z = 0; z < sizeZ; z += step) { // Increment z by step
                // Normalize the indices to the range [-1, 1]
                float posX = 2.0f * x / (sizeX - 1) - 1.0f;
                float posY = 2.0f * y / (sizeY - 1) - 1.0f;
                float posZ = 2.0f * z / (sizeZ - 1) - 1.0f;
                // Add position data
                this->vertexData_denoise.push_back(posX);
                this->vertexData_denoise.push_back(posY);
                this->vertexData_denoise.push_back(posZ);
                // Add color data (for example, based on Z value)
                QVector3D color = mapValueToColor(structure_dn[x][y][z], minVal, maxVal);
                this->vertexData_denoise.push_back(color.x()); // R
                this->vertexData_denoise.push_back(color.y()); // G
                this->vertexData_denoise.push_back(color.z()); // B
            }
        }
    }
    // The 2D QVector to store indices and values.
    QVector<QVector<double>> indices;
    QVector<QVector<double>> values;
    for (int i = 0; i < structure_dn.size(); ++i) {
        QVector<double> tempIndices;
        QVector<double> tempValues;
        for (int j = 0; j < structure_dn[i].size(); ++j) {
            for (int k = 0; k < structure_dn[i][j].size(); ++k) {
                // structure_dn[i][j][k] = structure_dn[i][j][k]/maxVal;
                if (structure_dn[i][j][k] > 0.7) {
                    tempIndices.append(static_cast<double>(k)); // Save the index along the 3rd dimension.
                    tempValues.append(structure_dn[i][j][k]); // Save the value.
                    break; // Move to the next 2D vector after finding the first value > 0.9.
                }
            }
            // If no value > 0.9 was found in the 3rd dimension, append -1 or another placeholder.
            if(tempIndices.size() == j) {
                tempIndices.append(-1); // Placeholder for not found.
                tempValues.append(0.0); // Placeholder value, assuming 0 is not a valid value in your context.
            }
        }
        indices.append(tempIndices);
        values.append(tempValues);
    }
    // // Debug output, replace with your method of saving or further processing.
    // qDebug() << "Indices:" << indices;
    // qDebug() << "Values:" << values;
    // Now save them to CSV files
    saveValuesToCSV(indices, "indices.csv");
    saveValuesToCSV(values, "values.csv");
    GraphWindow *graphWindow = new GraphWindow();
    graphWindow->setData(indices); // Set the data to be visualized
    graphWindow->show(); // Show the graph window
}

QVector3D MyOpenGLWidget::mapValueToColor(double value, double minVal, double maxVal) {
    // Normalize the value between 0 and 1
    double v = (value - minVal) / (maxVal - minVal);

    double r, g, b;
    r = g = b = 1.0; // Initialize with white color

    if (v < 0.25) {
        r = 0;
        g = 4 * v;
        b = 1;
    } else if (v < 0.5) {
        r = 0;
        g = 1;
        b = 1 - 4 * (v - 0.25);
    } else if (v < 0.75) {
        r = 4 * (v - 0.5);
        g = 1;
        b = 0;
    } else {
        r = 1;
        g = 1 - 4 * (v - 0.75);
        b = 0;
    }
    return QVector3D(r, g, b);
}

