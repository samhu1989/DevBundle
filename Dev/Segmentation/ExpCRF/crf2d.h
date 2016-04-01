#ifndef CRF2D_H
#define CRF2D_H
#include <QString>
#include <QObject>
#include <armadillo>
#include "common.h"
class CRF2D:public QObject
{
    Q_OBJECT
public:
    CRF2D(QImage& img,arma::uvec& lbl);
signals:
    int message(QString,int);
public slots:
    void process(void);
private:
    Config::Ptr config_;
    QImage& input_img_;
    arma::uvec& label_;
};

#endif // CRF2D_H
