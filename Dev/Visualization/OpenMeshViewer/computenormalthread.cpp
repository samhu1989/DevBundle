#include "computenormalthread.h"
#include "featurecore.h"
void ComputeNormalThread::run()
{
    Feature::computePointNormal(mesh_,0.0,int(21));
}

