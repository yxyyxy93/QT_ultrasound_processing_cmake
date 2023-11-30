#ifndef TRIMDIALOG_H
#define TRIMDIALOG_H

#include <QDialog>
#include <QSpinBox>

class TrimDialog : public QDialog {
    Q_OBJECT
public:
    explicit TrimDialog(QWidget *parent = nullptr);
    ~TrimDialog();

    int startI() const;
    int endI() const;
    int startJ() const;
    int endJ() const;
    int startK() const;
    int endK() const;

private:
    QSpinBox *startISpinBox;
    QSpinBox *endISpinBox;
    QSpinBox *startJSpinBox;
    QSpinBox *endJSpinBox;
    QSpinBox *startKSpinBox;
    QSpinBox *endKSpinBox;
};

#endif // TRIMDIALOG_H
