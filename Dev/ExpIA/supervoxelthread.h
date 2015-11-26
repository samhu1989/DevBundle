#ifndef SUPERVOXELTHREAD_H
#define SUPERVOXELTHREAD_H
#include <common.h>
#include <QThread>
#include "segmentationcore.h"
#include "voxelgraph.h"
class SupervoxelThread:public QThread
{
    Q_OBJECT
public:
    typedef Segmentation::SuperVoxelClustering<DefaultMesh> SvC;
    SupervoxelThread(MeshBundle<DefaultMesh>::PtrList&inputs):inputs_(inputs)
    {
        setObjectName("SupervoxelThread");
    }
signals:
    void message(QString,int);
public:
    bool configure(Config::Ptr config_);
protected:
    void run(void);
private:
    Segmentation::DefaultVoxelDistFunctor<DefaultMesh> vox_dist_;
    std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs_;
    Config::Ptr config_;

};

#endif // SUPERVOXELTHREAD_H
