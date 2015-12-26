#include "featureviewerwidget.h"
#include "ui_featureviewerwidget.h"
#include <QGLWidget>

FeatureViewerWidget::FeatureViewerWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::FeatureViewerWidget)
{
    ui->setupUi(this);
    ui->graphicsView->setViewport(new QGLWidget());
    ui->graphicsView->setScene(&scene_);
}

void FeatureViewerWidget::refresh()
{
    scene_.clear();
    feature_points_.clear();
    for(size_t idx;idx<features_.n_cols;++idx)
    {
        QPen pen;
        pen.setWidth(2);
        pen.setColor(QColor(qRgb(255,255,255)));
        QBrush brush;
        brush.setColor(
                    QColor(
                        qRgb(
                            feature_colors_(0,idx),
                            feature_colors_(1,idx),
                            feature_colors_(2,idx)
                            )
                        )
                    );
        QGraphicsEllipseItem* item = scene_.addEllipse(features_(0,idx),features_(1,idx),8.0,8.0,pen,brush);
        item->setToolTip(feature_strings_[idx]);
        feature_points_.push_back(item);
    }
}

FeatureViewerWidget::~FeatureViewerWidget()
{
    delete ui;
}
