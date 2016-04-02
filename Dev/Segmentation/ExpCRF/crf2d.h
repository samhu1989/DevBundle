#ifndef CRF2D_H
#define CRF2D_H
#include <QString>
#include <QObject>
#include <armadillo>
#include <QImage>
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
    QImage& input_img_;
    arma::uvec& label_;
};

#endif // CRF2D_H
