#ifndef MYOPENGLWIDGET_H
#define MYOPENGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLShaderProgram>
#include <QMatrix4x4>
#include <QMouseEvent>

class MyOpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT

public:
    explicit MyOpenGLWidget(QWidget *parent = nullptr);
    ~MyOpenGLWidget();
    void setData(const QVector<QVector<QVector<double>>> &data);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void prepareData(const QVector<QVector<QVector<double>>>& structure);
    QVector3D mapValueToColor(double value);

    void qNormalizeAngle(int &angle);
    void setXRotation(int angle);
    void setYRotation(int angle);
    void setZRotation(int angle);
    void updateViewMatrix();
    void mouseMoveEvent(QMouseEvent *event);

private:
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vbo;
    QOpenGLShaderProgram *m_program = nullptr;

    QVector<float> vertices; // Flattened array of vertex positions
    QVector<float> colors;   // Flattened array of vertex colors

    QMatrix4x4 projection;
    QMatrix4x4 view;

    int xRot = 0;
    int yRot = 0;
    int zRot = 0;
    QPoint lastPos; // Last position of the mouse


    // Shader source code would go here or be loaded from files
    const char *vertexShaderSource = R"(
        #version 330 core
        layout(location = 0) in vec3 aPos;
        layout(location = 1) in vec3 aColor;
        out vec3 ourColor;
        void main() {
            gl_Position = vec4(aPos, 1.0);
            ourColor = aColor;
        }
    )";

    const char *fragmentShaderSource = R"(
        #version 330 core
        in vec3 ourColor;
        out vec4 FragColor;
        void main() {
            FragColor = vec4(ourColor, 1.0);
        }
    )";
};

#endif // MYOPENGLWIDGET_H
