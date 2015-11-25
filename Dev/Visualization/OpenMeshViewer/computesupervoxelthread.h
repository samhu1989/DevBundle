#ifndef COMPUTESUPERVOXELTHREAD_H
#define COMPUTESUPERVOXELTHREAD_H
#include <QThread>
#include "common.h"
#include "segmentationcore.h"
class ComputeSupervoxelThread:public QThread
{
public:
    ComputeSupervoxelThread(MeshBundle<DefaultMesh>::Ptr in,MeshBundle<DefaultMesh>::Ptr out,QObject* parent=0):
        QThread(parent),
        input_(in),
        output_(out),
        svc_(0.05,0.19)
    {
        ;
    }
protected:
    void run();
private:
    MeshBundle<DefaultMesh>::Ptr input_;
    MeshBundle<DefaultMesh>::Ptr output_;
    Segmentation::SuperVoxelClustering<DefaultMesh> svc_;
};

#endif // COMPUTESUPERVOXELTHREAD_H
