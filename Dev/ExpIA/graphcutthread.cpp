#include "graphcutthread.h"
arma::sp_mat GraphCutThread::smooth_;

bool GraphCutThread::configure(Config::Ptr config)
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

void GraphCutThread::run(void)
{
    Segmentation::GraphCut gc;
    current_frame_ = 0;
    while( current_frame_ < meshes_.size() )
    {
        MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
        label_number_ = 1 + objects_.size();
        pix_number_ = m.graph_.voxel_centers.n_cols;
        gc.setLabelNumber( label_number_ );
        gc.setPixelNumber( pix_number_ );
        if(!prepareDataTerm())
        {
            std::cerr<<"Failed in prepareDataTerm"<<std::endl;
        }
        if(!current_data_||0==current_data_.use_count())std::logic_error("!current_data_||0==current_data_.use_count()");
        gc.inputDataTerm(current_data_);
        if(!prepareSmoothTerm(gc))
        {
            std::cerr<<"Failed in prepareSmoothTerm"<<std::endl;
        }
        gc.inputSmoothTerm(current_smooth_);
        gc.init(Segmentation::GraphCut::EXPANSION);
        if(!prepareNeighbors(gc))
        {
            std::cerr<<"Failed in prepareNeighbors"<<std::endl;
        }
        gc.updateInfo();
        float t;
        emit message(QString::fromStdString(gc.info()),0);
        gc.optimize(config_->getInt("GC_iter_num"),t);
        emit message(QString::fromStdString(gc.info()),0);
        arma::uvec sv_label;
        gc.getAnswer(sv_label);
        m.graph_.sv2pix(sv_label,outputs_[current_frame_]);
        m.custom_color_.fromlabel(outputs_[current_frame_]);
        ++current_frame_;
    }
}

void GraphCutThread::showMatch(size_t idx,DefaultMesh& mesh)
{
    MeshBundle<DefaultMesh>::Ptr m_ptr(new MeshBundle<DefaultMesh>);
    m_ptr->mesh_.request_vertex_colors();
    m_ptr->mesh_ = mesh;
    arma::Mat<uint8_t> cmat(
                (uint8_t*)m_ptr->custom_color_.vertex_colors(),
                4,
                m_ptr->mesh_.n_vertices(),
                false,
                true
    );
    cmat.fill(255);
    emit sendMatch(idx,m_ptr);
}

void GraphCutThread::showData(size_t current_label)
{
    if(current_label>=label_number_)return;
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    arma::mat data(data_.get(),label_number_,pix_number_,false,true);
    arma::Col<uint32_t> cv(pix_number_);
    arma::Col<uint32_t> cmat(
                (uint32_t*)m.custom_color_.vertex_colors(),
                m.mesh_.n_vertices(),
                false,
                true
    );
    float max_var = arma::max(data.row(current_label));
    float min_var = arma::min(data.row(current_label));
    float h;
    ColorArray::RGB32 tmp;
    int idx;
    #pragma omp for
    for(idx=0;idx<data.n_cols;++idx)
    {
        if( max_var!=min_var )h = ( data(current_label,idx) - min_var ) / ( max_var - min_var );
        else h = 0.0;
        ColorArray::hsv2rgb( h*255.0 + 5.0,0.5,1.0,tmp);
        cv(idx) = tmp.color;
    }
    m.graph_.sv2pix(cv,cmat);
}

void GraphCutThread::showSmooth()
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    arma::sp_mat r = arma::sum(smooth_);
    arma::sp_mat c = arma::sum(smooth_,1);
    arma::sp_mat var = c+r.t();
    arma::Col<uint32_t> cv(pix_number_);
    arma::Col<uint32_t> cmat(
                (uint32_t*)m.custom_color_.vertex_colors(),
                m.mesh_.n_vertices(),
                false,
                true
    );
    float max_var = arma::max(arma::max(var));
    float min_var = arma::min(arma::min(var));
    float h;
    ColorArray::RGB32 tmp;
    int idx;
    #pragma omp for
    for(idx=0;idx<var.size();++idx)
    {
        if( max_var!=min_var )h = ( var(idx) - min_var ) / ( max_var - min_var );
        else h = 0.0;
        ColorArray::hsv2rgb( h*255.0 + 5.0,0.5,1.0,tmp);
        cv(idx) = tmp.color;
    }
    m.graph_.sv2pix(cv,cmat);
}



bool GraphCutThread::prepareDataTerm()
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    data_.reset(new double[label_number_*pix_number_]);
    arma::mat data_mat(data_.get(),label_number_,pix_number_,false,true);
    std::vector<ObjModel::Ptr>::iterator oiter;
    uint32_t obj_index = 0;
    for(oiter=objects_.begin();oiter!=objects_.end();++oiter)
    {
        ObjModel& model = **oiter;
        DefaultMesh obj_mesh;
        if(model.transform(obj_mesh,current_frame_))
        {
            if(obj_mesh.n_vertices()==0)std::logic_error("obj_mesh.n_vertices()==0");
            if(config_->has("GC_show_match")&&1==config_->getInt("GC_show_match"))
            {
                showMatch(current_frame_,obj_mesh);
                QThread::sleep(1);
            }
            prepareDataForLabel(1+obj_index,m.graph_,obj_mesh,model.DistP_,model.NormP_,model.ColorP_);
            if(config_->has("GC_show_data")&&1==config_->getInt("GC_show_data"))
            {
                showData(1+obj_index);
                QThread::sleep(1);
            }
        }else{
            data_mat.row(1+obj_index).fill(0.0);
        }
        ++obj_index;
    }
    prepareDataForUnknown();
    if(config_->has("GC_show_data")&&1==config_->getInt("GC_show_data"))
    {
        showData(0);
        QThread::sleep(1);
    }
    normalizeData();
    current_data_.reset(new DataCost(data_.get()));
    return true;
}

void GraphCutThread::prepareDataForLabel(uint32_t l,
        VoxelGraph<DefaultMesh>& graph,
        DefaultMesh& obj,
        arma::fvec &dist_score,
        arma::fvec &norm_score,
        arma::fvec &color_score)
{
    double* data = data_.get();
    arma::mat data_mat(data,label_number_,pix_number_,false,true);
    arma::vec score;
    if(!config_->has("GC_color_var"))graph.match(obj,dist_score,norm_score,color_score,score,config_->getDouble("GC_distance_threshold"));
    else graph.match(obj,dist_score,norm_score,color_score,score,config_->getDouble("GC_distance_threshold"),config_->getDouble("GC_color_var"));
    data_mat.row(l) = score.t();
    if(!data_mat.row(l).is_finite())
    {
        std::cerr<<"infinite in data row "<<l<<std::endl;
    }
}

void GraphCutThread::prepareDataForUnknown()
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    arma::mat data((double*)data_.get(),label_number_,pix_number_,false,true);
    if(!data.is_finite())
    {
        std::cerr<<"infinite in data"<<std::endl;
    }
    arma::mat known_mat = data.rows(1,label_number_-1);
    arma::uword unknown_num = double(pix_number_)*config_->getDouble("GC_gamma_soft");
    std::cerr<<"unknown_num:"<<unknown_num<<std::endl;
    arma::rowvec known_median = arma::median(known_mat);
    arma::rowvec known_sum = 2.0*arma::sum(known_mat);
    arma::uvec sorted_i = arma::sort_index(known_sum);
    arma::rowvec unknown_data(data.n_cols);
    unknown_data = known_median;
    arma::uvec choosen_i = sorted_i.head(unknown_num);
    if(0!=choosen_i.size())unknown_data(choosen_i) = known_sum(choosen_i);
    data.row(0) = unknown_data;
}

void GraphCutThread::normalizeData()
{
    arma::mat data((double*)data_.get(),label_number_,pix_number_,false,true);
    arma::rowvec col_sum = arma::sum(data);
    #pragma omp for
    for(size_t i=0 ; i < pix_number_ ; ++i )
    {
        if( col_sum(i) !=0 )data.col(i) /= col_sum(i);
    }
    arma::uword unknown_num = double(pix_number_)*config_->getDouble("GC_gamma_hard");
    arma::rowvec unknown_data = data.row(0);
    arma::uvec sorted_i = arma::sort_index(unknown_data);
    arma::uvec choosen_i = sorted_i.head(unknown_num);
    unknown_data(choosen_i).fill(1.0);
    data = arma::mat(label_number_,pix_number_,arma::fill::ones) - data;
    data *= config_->getDouble("GC_data_weight");
    if(!data.is_finite())
    {
        std::cerr<<"infinite in data after normalized"<<std::endl;
    }
}

bool GraphCutThread::prepareSmoothTerm(Segmentation::GraphCut& gc)
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    smooth_ = arma::sp_mat(pix_number_,pix_number_);
    for( size_t idx = 0 ; idx < m.graph_.voxel_neighbors.n_cols ; ++idx )
    {
        size_t pix1 = m.graph_.voxel_neighbors(0,idx);
        size_t pix2 = m.graph_.voxel_neighbors(1,idx);
        if(pix1>=pix_number_)std::logic_error("pix1>=pix_number_");
        if(pix2>=pix_number_)std::logic_error("pix2>=pix_number_");
        if( pix1 < pix2 )
        {
            if(!config_->has("GC_color_var"))smooth_(pix1,pix2) = m.graph_.voxel_similarity(pix1,pix2);
            else smooth_(pix1,pix2) = m.graph_.voxel_similarity(pix1,pix2,config_->getDouble("GC_color_var"));
        }else{
            if(!config_->has("GC_color_var"))smooth_(pix2,pix1) = m.graph_.voxel_similarity(pix1,pix2);
            else smooth_(pix2,pix1) = m.graph_.voxel_similarity(pix1,pix2,config_->getDouble("GC_color_var"));
        }
    }
    smooth_ *= config_->getDouble("GC_smooth_weight");
    current_smooth_.reset(new SmoothnessCost(GraphCutThread::fnCost));
    if(config_->has("GC_show_smooth")&&1==config_->getInt("GC_show_smooth"))
    {
        showSmooth();
        if(config_->has("GC_show_smooth_sec"))
        {
            QThread::sleep(config_->getInt("GC_show_smooth_sec"));
        }else QThread::sleep(1);
    }
    return true;
}

MRF::CostVal GraphCutThread::fnCost(int pix1,int pix2,MRF::Label i,MRF::Label j)
{
    if(i==j)return 0.0;
    else if(pix1<pix2)return smooth_(pix1,pix2);
    else return smooth_(pix2,pix1);
}

bool GraphCutThread::prepareNeighbors(Segmentation::GraphCut& gc)
{
    MeshBundle<DefaultMesh>& m = *meshes_[current_frame_];
    double w_eps = 1.0 / double( m.mesh_.n_vertices() );
    for( size_t idx = 0 ; idx < m.graph_.voxel_neighbors.n_cols ; ++idx )
    {
        size_t pix1 = m.graph_.voxel_neighbors(0,idx);
        size_t pix2 = m.graph_.voxel_neighbors(1,idx);
        if(pix1>=pix_number_)std::logic_error("pix1>=pix_number_");
        if(pix2>=pix_number_)std::logic_error("pix2>=pix_number_");
        double w = w_eps * ( m.graph_.voxel_size(pix1) + m.graph_.voxel_size(pix2) );
        if(w<=0)std::logic_error("w<=0");
        if(!gc.setNeighbors(pix1,pix2,w))return false;
    }
    return true;
}
