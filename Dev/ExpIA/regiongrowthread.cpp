#include "regiongrowthread.h"
#include "segmentationcore.h"
#include "featurecore.h"
void RegionGrowThread::run()
{
    if(labels_.size()!=inputs_.size())
    {
        labels_.resize(inputs_.size());
    }
    InputIterator iiter;
    OutputIterator oiter;
    Segmentation::RegionGrowing<DefaultMesh> seg;

    seg.setNumberOfNeighbours(50);
    seg.setMinClusterSize(100);
    seg.setMaxClusterSize(std::numeric_limits<int>::max());
    seg.setSmoothnessThreshold(85.0/180.0*M_PI);
    seg.setCurvatureThreshold(std::numeric_limits<float>::max());

    oiter = labels_.begin();
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>& input = **iiter;
        emit message("Region Growing On: "+QString::fromStdString(input.name_),1000);
        std::shared_ptr<float> curvature;
        Feature::computePointNormal(input.mesh_,curvature,0.0,20);
        seg.setInputMesh(&input.mesh_);
        seg.getCurvatures() = curvature;
        seg.extract(*oiter);
        input.custom_color_.fromlabel(*oiter);
        ++oiter;
    }
}

