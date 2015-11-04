#ifndef COMPUTEOCTREETHREAD_H
#define COMPUTEOCTREETHREAD_H
#include <QThread>
#include "common.h"
class ComputeOctreeThread:public QThread
{
public:
    ComputeOctreeThread(DefaultMesh&in,DefaultMesh&out,QObject* parent=0):
        QThread(parent),
        input_(in),
        output_(out)
    {
        ;
    }
protected:
    void run();
private:
    DefaultMesh& input_;
    DefaultMesh& output_;
};

#endif // COMPUTEOCTREETHREAD_H
