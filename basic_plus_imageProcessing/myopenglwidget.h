#ifndef MYOPENGLWIDGET_H
#define MYOPENGLWIDGET_H

#include <QOpenGLWidget>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>

class MyOpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_3_Core {
    Q_OBJECT

public:
    explicit MyOpenGLWidget(QWidget *parent = nullptr, QString savename="");
    ~MyOpenGLWidget();
    void convertStructureToVertices(const QVector<QVector<QVector<double>>>& structure);

protected:
    void initializeGL() override;
    void resizeGL(int w, int h) override;
    void paintGL() override;
    void setViewAngle_x(int angle);
    void setViewAngle_y(int angle);
    void setVisibilityThreshold(int value);
    void toggleDenoisedData(int state);

    QVector3D mapValueToColor(double value, double minVal, double maxVal);

private:
    QOpenGLShaderProgram *m_program;
    GLuint VAO, VBO;

    QVector<GLfloat> vertexData;
    QVector<GLfloat> vertexData_denoise;

    double viewAngle_x;
    double viewAngle_y;
    double visibilityThreshold;
    bool useDenoisedData;

    QString fn;
};

#endif // MYOPENGLWIDGET_H
