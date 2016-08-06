#ifndef NCUT2D_H
#define NCUT2D_H
#include <QString>
#include <QObject>
#include <armadillo>
#include <QImage>
class NCut2D:public QObject
{
    Q_OBJECT
public:
    NCut2D(QImage& img,arma::uvec& lbl);
signals:
    int message(QString,int);
    void end();
public slots:
    void process(void);
private:
    QImage& input_img_;
    arma::uvec& label_;
};

#endif // NCUT2D_H
