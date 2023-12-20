#include "MyOpenGLWidget.h"
#include <QMouseEvent>

MyOpenGLWidget::MyOpenGLWidget(QWidget *parent)
    : QOpenGLWidget(parent), m_program(nullptr) {
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
    glBindVertexArray(0); // Unbind the VAO

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
    // Set the clear color to white
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // RGBA values for white
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPointSize(10.0f); // Set the point size to 10 pixels

    m_program->bind();

    // Calculate transformation
    QMatrix4x4 model;
    model.translate(0.0f, 0.0f, -2.0f); // Move back 2 units into the scene
    model.rotate(45.0f, 0.0f, 1.0f, 0.0f); // Rotate 45 degrees around the y-axis
    QMatrix4x4 view;
    view.translate(0.0f, 0.0f, -3.0f); // Move back 3 units from camera

    m_program->setUniformValue("model", model);
    m_program->setUniformValue("view", view);

    // Render the triangle
    glBindVertexArray(VAO);
    int numPoints = this->vertexData.size() / 6; // Assuming 3 position floats + 3 color floats per vertex
    glDrawArrays(GL_POINTS, 0, numPoints);

    // glBindVertexArray(0); // No need to unbind it every time
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

    // Now, minVal and maxVal can be used in the mapValueToColor function

    for (int x = 0; x < sizeX; ++x) {
        for (int y = 0; y < sizeY; ++y) {
            for (int z = 0; z < sizeZ; ++z) {
                // Normalize the indices to the range [-1, 1]
                float posX = 2.0f * x / (sizeX - 1) - 1.0f;
                float posY = 2.0f * y / (sizeY - 1) - 1.0f;
                float posZ = 2.0f * z / (sizeZ - 1) - 1.0f;
                // Add position data
                this->vertexData.push_back(posX);
                this->vertexData.push_back(posY);
                this->vertexData.push_back(posZ);
                // Add color data (for example, based on Z value)
                QVector3D color = mapValueToColor(structure[x][y][z], minVal, maxVal);
                this->vertexData.push_back(color.x()); // R
                this->vertexData.push_back(color.y()); // G
                this->vertexData.push_back(color.z()); // B
            }
        }
    }
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

