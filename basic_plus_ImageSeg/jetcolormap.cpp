#include "JetColorMap.h"
#include <QColor>

QColor JetColorMap::valueToJet(double normalizedValue) const {
    double r, g, b;

    if (normalizedValue < 0.125) {
        r = 0.0;
        g = 0.0;
        b = 4.0 * (normalizedValue + 0.125);
    } else if (normalizedValue < 0.375) {
        r = 0.0;
        g = 4.0 * (normalizedValue - 0.125);
        b = 1.0;
    } else if (normalizedValue < 0.625) {
        r = 4.0 * (normalizedValue - 0.375);
        g = 1.0;
        b = -4.0 * (normalizedValue - 0.625);
    } else if (normalizedValue < 0.875) {
        r = 1.0;
        g = -4.0 * (normalizedValue - 0.875);
        b = 0.0;
    } else {
        r = std::max(1.0, 4.0 * (normalizedValue - 0.875));
        g = 0.0;
        b = 0.0;
    }

    return QColor::fromRgbF(clamp(r, 0.0, 1.0), clamp(g, 0.0, 1.0), clamp(b, 0.0, 1.0));
}

double JetColorMap::clamp(double val, double min, double max) const {
    return std::max(min, std::min(max, val));
}
