#ifndef LOOPER_H
#define LOOPER_H
#include <QObject>
#include "common.h"
class Looper:public QObject
{
    Q_OBJECT
public:
    Looper(QObject* parent=0):QObject(parent){}
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
    void end();
public slots:
    void loop();
private:
    Config::Ptr config_;
};

#endif // LOOPER_H
