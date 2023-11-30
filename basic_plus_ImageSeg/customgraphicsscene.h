#ifndef CUSTOMGRAPHICSSCENE_H
#define CUSTOMGRAPHICSSCENE_H

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>

class CustomGraphicsScene : public QGraphicsScene
{
    Q_OBJECT
public:
    explicit CustomGraphicsScene(QObject *parent = nullptr);

    void mousePressEvent(QGraphicsSceneMouseEvent *event) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override;

    const QVector<QLineF>& getDrawnLines() const { return drawnLines; }
    QGraphicsLineItem* lineItemAt(const QPointF& p1, const QPointF& p2);
    void removeAllPolygons();
    void clearDrawnLines();

signals:
    void lineDrawn(const QLineF &line);

private:
    QGraphicsLineItem *tempLine = nullptr;
    bool drawing = false;
    QPointF startPoint;
    QPointF endPoint;

    QVector<QLineF> drawnLines;

};

#endif // CUSTOMGRAPHICSSCENE_H
