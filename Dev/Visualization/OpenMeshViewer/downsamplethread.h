#ifndef DOWNSAMPLETHREAD_H
#define DOWNSAMPLETHREAD_H
#include <QThread>
#include "common.h"
class DownSampleThread:public QThread
{
public:
    DownSampleThread(MeshBundle<DefaultMesh>::Ptr in,QObject* parent=0):
        QThread(parent),
        m_(in->mesh_)
    {
        ;
    }
protected:
    void run();
private:
    DefaultMesh& m_;
};

#endif // DOWNSAMPLETHREAD_H
