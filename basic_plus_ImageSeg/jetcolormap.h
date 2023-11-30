#ifndef JETCOLORMAP_H
#define JETCOLORMAP_H

#pragma once

#include <QColor>
#include <algorithm>

class JetColorMap {
public:
    QColor valueToJet(double normalizedValue) const;
private:
    double clamp(double val, double min, double max) const;
};

#endif // JETCOLORMAP_H
