#include "computesupervoxelthread.h"

using namespace Segmentation;
void ComputeSupervoxelThread::run()
{
    arma::uvec label;
    DefaultVoxelDistFunctor<DefaultMesh> distfunc;
    Segmentation::SuperVoxelClustering<DefaultMesh>::DistFunc dist;
    dist = std::bind(
                &DefaultVoxelDistFunctor<DefaultMesh>::dist,
                distfunc,
                std::placeholders::_1,
                std::placeholders::_2,
                svc_.getSeedResolution()
                );
    std::cerr<<"input:"<<std::endl;
    svc_.input(&input_->mesh_);
    std::cerr<<"set dist:"<<std::endl;
    svc_.setDistFunctor(dist);
    std::cerr<<"extract:"<<std::endl;
    svc_.extract(label);
    std::cerr<<"update color:"<<std::endl;
    input_->custom_color_.fromlabel(label);
}

