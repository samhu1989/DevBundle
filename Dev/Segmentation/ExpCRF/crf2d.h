#ifndef CRF2D_H
#define CRF2D_H
#include <QString>
#include <QObject>
class CRF2D:public QObject
{
    Q_OBJECT
public:
    CRF2D();
signals:
    int message(QString,int);
public slots:
    void process(void);
};

#endif // CRF2D_H
