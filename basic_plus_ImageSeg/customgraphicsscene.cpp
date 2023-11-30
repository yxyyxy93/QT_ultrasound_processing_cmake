#include "customgraphicsscene.h"
#include <QGraphicsLineItem>

CustomGraphicsScene::CustomGraphicsScene(QObject* parent) : QGraphicsScene(parent) {
    // Your constructor implementation
}

void CustomGraphicsScene::mousePressEvent(QGraphicsSceneMouseEvent *event) {
    this->startPoint = event->scenePos();
    tempLine = addLine(QLineF(this->startPoint, this->startPoint), QPen(Qt::green, 2));
}

void CustomGraphicsScene::mouseMoveEvent(QGraphicsSceneMouseEvent *event) {
    if (tempLine) {
        tempLine->setLine(QLineF(this->startPoint, event->scenePos()));
    }
}

void CustomGraphicsScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event) {
    this->endPoint = event->scenePos();
    if (tempLine) {
        tempLine->setLine(QLineF(startPoint, this->endPoint));
        emit lineDrawn(tempLine->line());
        tempLine = nullptr;
    }
    drawnLines.push_back(QLineF(startPoint, this->endPoint));
}

QGraphicsLineItem* CustomGraphicsScene::lineItemAt(const QPointF& p1, const QPointF& p2) {
    foreach (QGraphicsItem *item, items()) {
        if (QGraphicsLineItem *lineItem = qgraphicsitem_cast<QGraphicsLineItem *>(item)) {
            if (lineItem->line().p1() == p1 && lineItem->line().p2() == p2) {
                return lineItem;
            }
        }
    }
    return nullptr;
}

void CustomGraphicsScene::removeAllPolygons() {
    foreach (QGraphicsItem *item, items()) {
        if (QGraphicsPolygonItem *polygonItem = qgraphicsitem_cast<QGraphicsPolygonItem *>(item)) {
            removeItem(polygonItem);
            delete polygonItem;
        }
    }
}


void CustomGraphicsScene::clearDrawnLines() {
    drawnLines.clear();
}
