#include "supervoxelthread.h"
bool SupervoxelThread::configure(Config::Ptr config)
{
    config_ = config;
    if(inputs_.empty())return false;
    if(!config_->has("Sv_resolution"))return false;
    if(!config_->has("Sv_seed_resolution"))return false;
    if(!config_->has("Sv_spatial_weight"))return false;
    if(!config_->has("Sv_color_weight"))return false;
    if(!config_->has("Sv_normal_weight"))return false;
    return true;
}

void SupervoxelThread::run(void)
{
    SvC svc(
            config_->getFloat("Sv_seed_resolution"),
            config_->getFloat("Sv_resolution")
        );
    SvC::DistFunc dist = std::bind(
                &Segmentation::DefaultVoxelDistFunctor<DefaultMesh>::dist,
                vox_dist_,
                std::placeholders::_1,
                std::placeholders::_2,
                svc.getSeedResolution()
                );
    svc.setDistFunctor(dist);
    timer_.restart();
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>::Ptr &input = *iter;
        svc.input(&input->mesh_);
        svc.extract(input->graph_.voxel_label);
        input->custom_color_.fromlabel(input->graph_.voxel_label);
        svc.getCentroids(input->graph_.voxel_centers);
        svc.getCentroidColors(input->graph_.voxel_colors);
        svc.getCentroidNormals(input->graph_.voxel_normals);
        svc.getSizes(input->graph_.voxel_size);
        svc.getSupervoxelAdjacency(input->graph_.voxel_neighbors);
    }
    QString msg;
    int ms = timer_.elapsed();
    int s = ms/1000;
    ms -= s*1000;
    int m = s/60;
    s -= m*60;
    int h = m/60;
    m -= h*60;
    msg = msg.sprintf("Time Used of Region Grow:%2u:%2u:%2u.%3u",h,m,s,ms);
    emit message(msg,0);
}


