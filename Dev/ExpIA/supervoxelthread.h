#ifndef SUPERVOXELTHREAD_H
#define SUPERVOXELTHREAD_H
#include <QThread>
class SupervoxelThread:public QThread
{
    Q_OBJECT
public:
    SupervoxelThread();
signals:
    void message(QString,int);
public slots:

protected:
    void run(void);

};

#endif // SUPERVOXELTHREAD_H
