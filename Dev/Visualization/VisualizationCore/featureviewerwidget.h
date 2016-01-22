#ifndef FEATUREVIEWERWIDGET_H
#define FEATUREVIEWERWIDGET_H
#include "visualizationcore_global.h"
#include <QWidget>
#include <featureviewerwidget.h>
#include "common.h"
#include <QGraphicsScene>
#include <QGraphicsEllipseItem>
#include <QWheelEvent>
namespace Ui {
class FeatureViewerWidget;
}

class VISUALIZATIONCORESHARED_EXPORT FeatureViewerWidget : public QWidget
{
    Q_OBJECT

public:
    explicit FeatureViewerWidget(QWidget *parent = 0);
    ~FeatureViewerWidget();
    inline void set_features(const arma::fmat& f){features_ = f;}
    inline void set_feature_colors(const arma::Mat<uint8_t>& color){feature_colors_ = color;}
    inline void set_feature_string(const QStringList& string){feature_strings_ = string;}
public slots:
    virtual void refresh();
protected:
    void wheelEvent(QWheelEvent*);
protected:
    arma::fmat features_;
    arma::Mat<uint8_t> feature_colors_;
    QStringList feature_strings_;
    QList<QGraphicsEllipseItem*> feature_points_;
private:
    Ui::FeatureViewerWidget *ui;
    QGraphicsScene scene_;
    qreal scale_;
};

#endif // FEATUREVIEWERWIDGET_H