#ifndef COMPUTESUPERVOXELTHREAD_H
#define COMPUTESUPERVOXELTHREAD_H
#include <QThread>
#include "common.h"
class ComputeSupervoxelThread:public QThread
{
public:
    ComputeSupervoxelThread(DefaultMesh&in,DefaultMesh&out,QObject* parent=0):
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

#endif // COMPUTESUPERVOXELTHREAD_H
