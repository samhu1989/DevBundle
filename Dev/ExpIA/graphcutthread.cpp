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
        if(!prepareDataTerm())
        {

        }
        if(!prepareSmoothTerm())
        {

        }
        gc.setLabelNumber( 1 + objects_.size() );
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
    ;
}

bool GraphCutThread::prepareSmoothTerm()
{
    ;
}
