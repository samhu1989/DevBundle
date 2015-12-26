#include "inpatchgraphcut.h"

bool InPatchGraphCut::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("GC_iter_num"))return false;
    if(!config_->has("GC_data_weight"))return false;
    if(!config_->has("GC_smooth_weight"))return false;
    if(!config_->has("GC_gamma_soft"))return false;
    if(config_->getDouble("GC_gamma_soft")>0.5||config_->getDouble("GC_gamma_soft")<0){
        emit message(tr("Wrong GC_gamma_soft"),0);
        return false;
    }
    if(!config_->has("GC_gamma_hard"))return false;
    if(config_->getDouble("GC_gamma_hard")>0.5||config_->getDouble("GC_gamma_hard")<0){
        emit message(tr("Wrong GC_gamma_hard"),0);
        return false;
    }
    if(!config_->has("GC_distance_threshold"))return false;
    if(meshes_.empty())return false;
    if(meshes_.front()->graph_.empty())return false;
    if(objects_.empty())return false;
    return true;
}

void InPatchGraphCut::run(void)
{
    current_frame_ = 0;
    current_patch_mesh_.request_vertex_colors();
    current_patch_mesh_.request_vertex_normals();
    while( current_frame_ < meshes_.size() )
    {
        timer.restart();
        std::cerr<<"Frame "<<current_frame_<<std::endl;
        for_each_frame();
        QString msg;
        msg = msg.sprintf("%u ms for F %u",timer.elapsed(),current_frame_);
        emit message(msg,0);
        ++current_frame_;
    }
}

void InPatchGraphCut::for_each_frame(void)
{
    Segmentation::GraphCut gc;
    MeshBundle<DefaultMesh>& mesh = *meshes_[current_frame_];
    arma::uvec& label = labels_[current_frame_];
    arma::uword label_max = arma::max(label);
    current_label_ = 1;
    while( current_label_ <= label_max )
    {
        arma::uvec label_indices = arma::find( label == current_label_ );
        if(label_indices.is_empty()){
            ++ current_label_;
            continue;
        }
//        std::cerr<<"extract patch graph"<<std::endl;
        current_patch_mesh_.clear();
        current_patch_graph_ = VoxelGraph<DefaultMesh>::getSubGraphPtr(mesh.graph_,label_indices,current_patch_mesh_);
//        std::cerr<<"done extract patch graph"<<std::endl;
        label_number_ = 2;
        pix_number_ = current_patch_graph_->voxel_centers.n_cols;
        gc.setLabelNumber( label_number_ );
        gc.setPixelNumber( pix_number_ );
//        std::cerr<<"prepareDataTerm()"<<std::endl;
        if(!prepareDataTerm(gc))
        {
            std::cerr<<"Failed to prepare data term"<<std::endl;
        }
//        std::cerr<<"prepareSmoothTerm()"<<std::endl;
        if(!prepareSmoothTerm(gc))
        {
            std::cerr<<"Failed to prepare smooth term"<<std::endl;
        }
        gc.init(Segmentation::GraphCut::EXPANSION);
//        std::cerr<<"prepareNeighbors()"<<std::endl;
        if(!prepareNeighbors(gc))
        {
            std::cerr<<"Failed to prepare graph"<<std::endl;
        }
//        std::cerr<<"gc.updateInfo()"<<std::endl;
        gc.updateInfo();
        float t;
        emit message(QString::fromStdString(gc.info()),0);
        gc.optimize(config_->getInt("GC_iter_num"),t);
        emit message(QString::fromStdString(gc.info()),0);
        arma::uvec gc_label;
        gc.getAnswer(gc_label);
//        std::cerr<<"get answer"<<std::endl;
        applyToFrame(gc_label,label_indices);
//        std::cerr<<"done apply to frame"<<std::endl;
        ++current_label_;
    }
}

bool InPatchGraphCut::prepareDataTerm(Segmentation::GraphCut& gc)
{
//    std::cerr<<"0"<<std::endl;
    data_.reset(new double[label_number_*pix_number_]);
    double* data = data_.get();
    arma::mat data_mat(data,label_number_,pix_number_,false,true);
    ObjModel& model = *objects_[ current_label_ - 1 ];
    DefaultMesh obj_mesh;
    model.transform(obj_mesh,current_frame_);
    arma::vec score;
//    std::cerr<<"1"<<std::endl;
    if(!config_->has("GC_color_var"))current_patch_graph_->match(obj_mesh,model.DistP_,model.NormP_,model.ColorP_,score,config_->getDouble("GC_distance_threshold"));
    else current_patch_graph_->match(obj_mesh,model.DistP_,model.NormP_,model.ColorP_,score,config_->getDouble("GC_distance_threshold"),config_->getDouble("GC_color_var"));
//    std::cerr<<"2"<<std::endl;
    data_mat.row(0) = score.t();
    data_mat.row(1).fill(1.0);
    data_mat.row(1) -= data_mat.row(0);
    if(!data_mat.is_finite())
    {
        std::cerr<<"infinite in data of frame "<<current_frame_<<" patch "<<current_label_<<std::endl;
    }
//    std::cerr<<"3"<<std::endl;
    current_data_.reset(new DataCost(data_.get()));
    gc.inputDataTerm(current_data_);
    return true;
}

arma::sp_mat InPatchGraphCut::smooth_;
MRF::CostVal InPatchGraphCut::fnCost(int pix1,int pix2,MRF::Label i,MRF::Label j)
{
    if(i==j)return 0.0;
    else if(pix1<pix2)return smooth_(pix1,pix2);
    else return smooth_(pix2,pix1);
}

bool InPatchGraphCut::prepareSmoothTerm(Segmentation::GraphCut& gc)
{
    smooth_ = arma::sp_mat(pix_number_,pix_number_);
    for( size_t idx = 0 ; idx < current_patch_graph_->voxel_neighbors.n_cols ; ++idx )
    {
        size_t pix1 = current_patch_graph_->voxel_neighbors(0,idx);
        size_t pix2 = current_patch_graph_->voxel_neighbors(1,idx);
        if(pix1>=pix_number_)std::logic_error("pix1>=pix_number_");
        if(pix2>=pix_number_)std::logic_error("pix2>=pix_number_");
        if( pix1 < pix2 )
        {
            if(!config_->has("GC_color_var"))smooth_(pix1,pix2) = current_patch_graph_->voxel_similarity(pix1,pix2);
            else smooth_(pix1,pix2) = current_patch_graph_->voxel_similarity(pix1,pix2,config_->getDouble("GC_color_var"));
        }else{
            if(!config_->has("GC_color_var"))smooth_(pix2,pix1) = current_patch_graph_->voxel_similarity(pix1,pix2);
            else smooth_(pix2,pix1) = current_patch_graph_->voxel_similarity(pix1,pix2,config_->getDouble("GC_color_var"));
        }
    }
    smooth_ *= config_->getDouble("GC_smooth_weight");
    current_smooth_.reset(new SmoothnessCost(InPatchGraphCut::fnCost));
    gc.inputSmoothTerm(current_smooth_);
    return true;
}

bool InPatchGraphCut::prepareNeighbors(Segmentation::GraphCut& gc)
{
    double w_eps = 1.0 / double( current_patch_mesh_.n_vertices() );
    for( size_t idx = 0 ; idx < current_patch_graph_->voxel_neighbors.n_cols ; ++idx )
    {
        size_t pix1 = current_patch_graph_->voxel_neighbors(0,idx);
        size_t pix2 = current_patch_graph_->voxel_neighbors(1,idx);
        if(pix1>=pix_number_)std::logic_error("pix1>=pix_number_");
        if(pix2>=pix_number_)std::logic_error("pix2>=pix_number_");
        double w = w_eps * ( current_patch_graph_->voxel_size(pix1) + current_patch_graph_->voxel_size(pix2) );
        if(w<=0)std::logic_error("w<=0");
        if(!gc.setNeighbors(pix1,pix2,w))return false;
    }
    return true;
}

void InPatchGraphCut::applyToFrame(const arma::uvec& gc_label,const arma::uvec& label_indices)
{
    MeshBundle<DefaultMesh>& mesh = *meshes_[current_frame_];
    arma::uvec& label = labels_[current_frame_];
    arma::uvec patch_label;
    current_patch_graph_->sv2pix(gc_label,patch_label);
    if(patch_label.size()!=label_indices.size())std::logic_error("gc_label.size()!=label_indices.size()");
    arma::uvec new_patch_indices = arma::find( patch_label==0 );
    arma::uvec new_patch_value = label(label_indices);
    new_patch_value(new_patch_indices).fill(0);
    label(label_indices) = new_patch_value;
    mesh.custom_color_.fromlabel(label);
}

