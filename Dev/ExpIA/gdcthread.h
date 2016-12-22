#ifndef GDCTHREAD_H
#define GDCTHREAD_H
#include "common.h"
#include <QObject>
#include <QTimer>
class GDCThread:public QObject
{
    Q_OBJECT
public:
    explicit GDCThread(
            std::vector<QWidget*> &inputs,
            QTimer &gl_timer,
            QObject* parent=0
            );
    bool configure(Config::Ptr);
signals:
    void message(QString,int);
    void end();
public slots:
    void process();
protected:
    void feature_to_color(uint32_t* ptr,arma::uword size,const arma::fmat& f);
private:
    std::vector<QWidget*>& inputs_;
    Config::Ptr config_;
};

#endif // GDCTHREAD_H
