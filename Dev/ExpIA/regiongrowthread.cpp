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

    seg.setNumberOfNeighbours(config_->getInt("RegionGrow_k"));
    seg.setMinClusterSize(config_->getInt("RegionGrow_cluster_min_num"));
    seg.setMaxClusterSize(config_->getInt("RegionGrow_cluster_max_num"));
    seg.setSmoothnessThreshold(config_->getFloat("RegionGrow_max_norm_angle")/180.0*M_PI);
    seg.setCurvatureThreshold(std::numeric_limits<float>::max());

    oiter = labels_.begin();
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>& input = **iiter;
        emit message("Region Growing On: "+QString::fromStdString(input.name_),1000);
        std::shared_ptr<float> curvature;
        Feature::computePointNormal(input.mesh_,curvature,0.0,config_->getInt("NormalEstimation_k"));
        seg.setInputMesh(&input.mesh_);
        arma::uvec indices = arma::find(*oiter==0);
        seg.setIndices(indices);
        seg.getCurvatures() = curvature;
        seg.extract(*oiter);
        input.custom_color_.fromlabel(*oiter);
        ++oiter;
    }
}

bool RegionGrowThread::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("NormalEstimation_k"))return false;
    if(!config_->has("RegionGrow_k"))return false;
    if(!config_->has("RegionGrow_cluster_min_num"))return false;
    if(!config_->has("RegionGrow_cluster_max_num"))return false;
    if(!config_->has("RegionGrow_max_norm_angle"))return false;
    return true;
}

