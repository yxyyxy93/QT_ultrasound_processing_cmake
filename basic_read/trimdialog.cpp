#include "trimdialog.h"

#include <QDialog>
#include <QSpinBox>
#include <QFormLayout>
#include <QPushButton>
#include <QDialogButtonBox>

TrimDialog::TrimDialog(QWidget *parent) : QDialog(parent) {
    QFormLayout *layout = new QFormLayout(this);
    // Create spin boxes for the start and end indices
    startISpinBox = new QSpinBox(this);
    endISpinBox = new QSpinBox(this);
    startJSpinBox = new QSpinBox(this);
    endJSpinBox = new QSpinBox(this);
    startKSpinBox = new QSpinBox(this);
    endKSpinBox = new QSpinBox(this);

    // Set some reasonable maximum values for the spin boxes
    startISpinBox->setMaximum(1000);
    endISpinBox->setMaximum(1000);
    startJSpinBox->setMaximum(1000);
    endJSpinBox->setMaximum(1000);
    startKSpinBox->setMaximum(1000);
    endKSpinBox->setMaximum(10000);

    // Add the spin boxes to the layout
    layout->addRow("Start I", startISpinBox);
    layout->addRow("End I", endISpinBox);
    layout->addRow("Start J", startJSpinBox);
    layout->addRow("End J", endJSpinBox);
    layout->addRow("Start K", startKSpinBox);
    layout->addRow("End K", endKSpinBox);

    // Add OK and Cancel buttons
    QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    connect(buttonBox, &QDialogButtonBox::accepted, this, &TrimDialog::accept);
    connect(buttonBox, &QDialogButtonBox::rejected, this, &TrimDialog::reject);
    layout->addWidget(buttonBox);
}

TrimDialog::~TrimDialog() {
    // Any necessary cleanup here...
}


int TrimDialog::startI() const { return startISpinBox->value(); }
int TrimDialog::endI() const { return endISpinBox->value(); }
int TrimDialog::startJ() const { return startJSpinBox->value(); }
int TrimDialog::endJ() const { return endJSpinBox->value(); }
int TrimDialog::startK() const { return startKSpinBox->value(); }
int TrimDialog::endK() const { return endKSpinBox->value(); }
