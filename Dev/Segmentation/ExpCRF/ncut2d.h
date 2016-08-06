#ifndef NCUT2D_H
#define NCUT2D_H
#include <QString>
#include <QObject>
#include <armadillo>
#include <QImage>
#include "common.h"
#include "segmentationcore.h"
class NCut2D:public QObject
{
    Q_OBJECT
public:
    NCut2D(QImage& img,arma::uvec& lbl);
    bool configure(Config::Ptr);
signals:
    int message(QString,int);
    void end();
public slots:
    void process(void);
private:
    Segmentation::NormalizedCuts<DefaultMesh> cut_;
    QImage& input_img_;
    arma::uvec& label_;
};

#endif // NCUT2D_H
