#ifndef MYOPENGLWIDGET_H
#define MYOPENGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class MyOpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT

public:
    explicit MyOpenGLWidget(QWidget *parent = nullptr);
    ~MyOpenGLWidget();
    void convertStructureToVertices(const QVector<QVector<QVector<double>>>& structure);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;

    QVector3D mapValueToColor(double value, double minVal, double maxVal);

private:
    QOpenGLShaderProgram *m_program;
    GLuint VAO, VBO;

    QVector<GLfloat> vertexData;
};

#endif // MYOPENGLWIDGET_H
