#include "MyOpenGLWidget.h"

MyOpenGLWidget::MyOpenGLWidget(QWidget *parent)
    : QOpenGLWidget(parent) {
}

MyOpenGLWidget::~MyOpenGLWidget() {
    makeCurrent();
    m_vbo.destroy();
    m_vao.destroy();
    delete m_program;
    doneCurrent();
}

void MyOpenGLWidget::initializeGL() {
    initializeOpenGLFunctions();
    m_program = new QOpenGLShaderProgram(this);
    // This should be removed as it's duplicating the work done below.
    m_program->addShaderFromSourceCode(QOpenGLShader::Vertex, vertexShaderSource);
    m_program->addShaderFromSourceCode(QOpenGLShader::Fragment, fragmentShaderSource);
    m_program->link();
    m_program->bind();
    // Link shader program
    if (!m_program->link()) {
        qWarning() << "Shader Link Error:" << m_program->log();
        return; // Don't proceed further if there's a linking error
    }
    // Define your vertex data
    GLfloat vertices[] = {
        // positions         // colors
        0.5f, -0.5f, 0.0f,  1.0f, 0.0f, 0.0f,
        -0.5f, -0.5f, 0.0f,  0.0f, 1.0f, 0.0f,
        0.0f,  0.5f, 0.0f,  0.0f, 0.0f, 1.0f
    };
    m_vao.create();
    m_vao.bind();
    m_vbo.create();
    m_vbo.bind();
    m_vbo.allocate(vertices, sizeof(vertices));
    // Position attribute
    m_program->enableAttributeArray(0);
    m_program->setAttributeBuffer(0, GL_FLOAT, 0, 3, sizeof(GLfloat) * 6);
    // Color attribute
    m_program->enableAttributeArray(1);
    m_program->setAttributeBuffer(1, GL_FLOAT, sizeof(GLfloat) * 3, 3, sizeof(GLfloat) * 6);
    m_vao.release();
    m_vbo.release();
    m_program->release();
    // Set up the view matrix for the camera
    // Adjust these values as needed for your specific data and desired view
    QVector3D cameraPos(2.0f, 2.0f, 2.0f); // Position the camera in 3D space for a perspective view
    QVector3D cameraFocus(0.0f, 0.0f, 0.0f); // Look at the center of the data
    QVector3D cameraUp(0.0f, 0.0f, 1.0f); // Z is up in this coordinate system
    view.lookAt(cameraPos, cameraFocus, cameraUp);
    glEnable(GL_DEPTH_TEST);
}

void MyOpenGLWidget::resizeGL(int w, int h) {
    // Update the projection matrix
    projection.setToIdentity();
    float aspectRatio = float(w) / float(h);
    projection.perspective(45.0f, aspectRatio, 0.1f, 10.0f);
}

void MyOpenGLWidget::paintGL() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    m_program->bind();
    m_program->setUniformValue("projection", projection);
    m_program->setUniformValue("view", view);
    m_vao.bind();
    // Set the point size
    glPointSize(5.0f); // Adjust the size as needed
    // Assuming you've already populated the vertices array
    m_vbo.bind();
    m_vbo.allocate(vertices.constData(), vertices.size() * sizeof(float));
    glDrawArrays(GL_POINTS, 0, vertices.size() / 6); // 6 = 3 for position + 3 for color
    m_vao.release();
    m_program->release();
}

void MyOpenGLWidget::prepareData(const QVector<QVector<QVector<double>>>& structure) {
    // Clear any previous data
    vertices.clear();
    // Find min and max values for normalization
    double minVal = std::numeric_limits<double>::max();
    double maxVal = std::numeric_limits<double>::lowest();
    for (const auto& layer : structure) {
        for (const auto& row : layer) {
            for (double val : row) {
                minVal = std::min(minVal, val);
                maxVal = std::max(maxVal, val);
            }
        }
    }
    const int sizeX = structure.size();
    const int sizeY = structure.size() > 0 ? structure[0].size() : 0;
    const int sizeZ = structure.size() > 0 && structure[0].size() > 0 ? structure[0][0].size() : 0;
    for (int x = 0; x < sizeX; ++x) {
        for (int y = 0; y < sizeY; ++y) {
            for (int z = 0; z < sizeZ; ++z) {
                // Normalize the x, y, z coordinates
                float normalizedX = (static_cast<float>(x) / (sizeX - 1)) * 2.0f - 1.0f;
                float normalizedY = (static_cast<float>(y) / (sizeY - 1)) * 2.0f - 1.0f;
                float normalizedZ = (static_cast<float>(z) / (sizeZ - 1)) * 2.0f - 1.0f;
                // Normalize the value
                double normalizedValue = (structure[x][y][z] - minVal) / (maxVal - minVal);
                QVector3D color = mapValueToColor(normalizedValue);
                // Position
                vertices.push_back(normalizedX);
                vertices.push_back(normalizedY);
                vertices.push_back(normalizedZ);
                // Color
                vertices.push_back(color.x());
                vertices.push_back(color.y());
                vertices.push_back(color.z());
            }
        }
    }
}

void MyOpenGLWidget::setData(const QVector<QVector<QVector<double>>> &data) {
    prepareData(data); // Call the existing method to prepare the data
    update();          // Trigger a repaint
}

QVector3D MyOpenGLWidget::mapValueToColor(double value) {
    double r, g, b;
    // This is a simple implementation of the Jet colormap.
    const double v = std::max(0.0, std::min(value, 1.0));

    if (v <= 0.25) {
        r = 0.0;
        g = 4 * v;
        b = 1.0;
    } else if (v <= 0.5) {
        r = 0.0;
        g = 1.0;
        b = 1.0 + 4 * (0.25 - v);
    } else if (v <= 0.75) {
        r = 4 * (v - 0.5);
        g = 1.0;
        b = 0.0;
    } else {
        r = 1.0;
        g = 1.0 - 4 * (v - 0.75);
        b = 0.0;
    }

    return QVector3D(r, g, b);
}


// ****************** mouse interaction
// Helper function to normalize angles
void MyOpenGLWidget::qNormalizeAngle(int &angle) {
    while (angle < 0)
        angle += 360 * 16;
    while (angle > 360 * 16)
        angle -= 360 * 16;
}

// Setters for rotation
void MyOpenGLWidget::setXRotation(int angle) {
    qNormalizeAngle(angle);
    if (angle != xRot) {
        xRot = angle;
        updateViewMatrix();
        update();
    }
}

void MyOpenGLWidget::setYRotation(int angle) {
    qNormalizeAngle(angle);
    if (angle != yRot) {
        yRot = angle;
        updateViewMatrix();
        update();
    }
}

void MyOpenGLWidget::setZRotation(int angle) {
    qNormalizeAngle(angle);
    if (angle != zRot) {
        zRot = angle;
        updateViewMatrix();
        update();
    }
}

// This method updates the view matrix based on the current rotation angles.
void MyOpenGLWidget::updateViewMatrix() {
    QMatrix4x4 rotation;
    rotation.rotate(xRot / 16.0f, 1.0f, 0.0f, 0.0f);
    rotation.rotate(yRot / 16.0f, 0.0f, 1.0f, 0.0f);
    rotation.rotate(zRot / 16.0f, 0.0f, 0.0f, 1.0f);
    // Assuming your initial view matrix is set up with lookAt.
    QVector3D cameraPos(3.0f, 3.0f, 3.0f);
    QVector3D cameraUp(0.0f, 0.0f, 1.0f);
    QVector3D cameraFocus(0.0f, 0.0f, 0.0f);
    view.setToIdentity();
    view.lookAt(cameraPos, cameraFocus, cameraUp);
    view = view * rotation; // Apply the rotation to the view matrix
}

// call the setXRotation, setYRotation, and setZRotation methods
void MyOpenGLWidget::mouseMoveEvent(QMouseEvent *event) {
    int dx = event->x() - lastPos.x();
    int dy = event->y() - lastPos.y();
    if (event->buttons() & Qt::LeftButton) {
        // Adjust these values if the rotation is too fast or too slow
        setXRotation(xRot + 8 * dy);
        setYRotation(yRot + 8 * dx);
    } else if (event->buttons() & Qt::RightButton) {
        // Right mouse button drag could be used for z-axis rotation or other actions like zoom
        setZRotation(zRot + 8 * dx);
    }
    lastPos = event->pos(); // Update the position for the next event
}
