#include "graphcutthread.h"
bool GraphCutThread::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("GC_n_iter"))return false;
    return true;
}

void GraphCutThread::run(void)
{
    Segmentation::GraphCut gc;
    current_frame_ = 0;
    while( current_frame_ < meshes_.size() )
    {
        MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
        label_number_ = 1 + objects_.size();
        pix_number_ = m.graph_.voxel_centers.n_cols;
        if(!prepareDataTerm())
        {

        }
        if(!prepareSmoothTerm())
        {

        }
        gc.setLabelNumber( label_number_ );
        gc.inputDataTerm(current_data_);
        gc.inputSmoothTerm(current_smooth_);
        gc.init(m.graph_,Segmentation::GraphCut::EXPANSION);
        float t;
        gc.optimize(config_->getInt("GC_n_iter"),t);
        emit message(QString::fromStdString(gc.info()),0);
        arma::uvec sv_label;
        gc.getAnswer(sv_label);
        m.graph_.sv2pix(sv_label,outputs_[current_frame_]);
        ++current_frame_;
    }
}

bool GraphCutThread::prepareDataTerm()
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    data_.reset(new double[label_number_*pix_number_]);
    std::vector<ObjModel::Ptr>::iterator oiter;
    uint32_t obj_index = 0;
    for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
    {
        ObjModel& model = **oiter;
        DefaultMesh obj_mesh;
        model.transform(obj_mesh,current_frame_);
        prepareDataForLabel(1+obj_index,m.graph_,obj_mesh,model.GeoP_,model.ColorP_);
        ++obj_index;
    }
    current_data_.reset(new DataCost(data_.get()));
}

void GraphCutThread::prepareDataForLabel(uint32_t l,
        VoxelGraph& graph,
        DefaultMesh& obj,
        std::vector<float> &geo_score,
        std::vector<float> &color_score)
{
    double* data = data_.get();
    std::vector<float> score;
    graph.match(obj,geo_score,color_score,score);
    for(size_t pix=0 ; pix < pix_number_ ; ++ pix)
    {
        data[ pix*label_number_ + l ] = score[pix];
    }
}

bool GraphCutThread::prepareSmoothTerm()
{
    ;
}
