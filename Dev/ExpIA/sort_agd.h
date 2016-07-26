#ifndef SORT_AGD_H
#define SORT_AGD_H
#include <QObject>
class Sort_AGD:public QObject
{
    Q_OBJECT
public:
    Sort_AGD();
public slots:
    void process();
protected:
private:
};

#endif // SORT_AGD_H
