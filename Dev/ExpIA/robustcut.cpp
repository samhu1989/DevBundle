#include "robustcut.h"
#include <QMessageBox>
std::vector<arma::umat> RobustCut::base_segment_list_;
RobustCut::RobustCut(
        MeshBundle<DefaultMesh>::PtrList &inputs,
        std::vector<arma::uvec> &labels,
        QObject *parent
        ):inputs_(inputs),labels_(labels),QObject(parent)
{
    setObjectName("RobustCut");
    if(!base_segment_list_.empty())
    {
        base_segment_i_ = 0;
    }
    base_segment_N_ = std::numeric_limits<uint32_t>::max();
    for(std::vector<arma::umat>::iterator iter=base_segment_list_.begin();iter!=base_segment_list_.end();++iter)
    {
        arma::uword n = iter->n_cols;
        if(n<base_segment_N_)base_segment_N_=n;
    }
}

bool RobustCut::configure(Config::Ptr config)
{
    if(inputs_.empty())return false;
    if(inputs_.front()->graph_.empty()){
        QMessageBox::warning(NULL,tr("RobustCut"),tr("Please build graph by supervoxel first!"));
        return false;
    }
    return cuts_.configure(config);
}

void RobustCut::base_segments()
{
    base_segment_list_.clear();
    base_segment_i_ = 0;
    base_segment_N_ = std::numeric_limits<uint32_t>::max();
    std::vector<arma::uvec>::iterator oiter;
    oiter = labels_.begin();
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    arma::uword i = 1;
    QString msg;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        std::cerr<<"============"<<std::endl;
        cuts_.generate_base_segments(mesh_ptr,*oiter,true);
        std::cerr<<"============"<<std::endl;
        base_segment_list_.emplace_back(cuts_.base_segments());
        arma::uword n = base_segment_list_.back().n_cols;
        if(n<base_segment_N_)base_segment_N_=n;
        ++oiter;
        if(oiter==labels_.end())break;
        emit message(msg.sprintf("Done %u/%u",i,inputs_.size()),0);
        ++ i;
    }
    show_base_segment();
}

void RobustCut::keyPressEvent(QKeyEvent* event)
{
    if(event->type()!=QEvent::KeyPress)return;
    if(event->key()==Qt::Key_Escape)exit_base_segment();
    if(event->key()==Qt::Key_Left)last_base_segment();
    if(event->key()==Qt::Key_Right)next_base_segment();
}

void RobustCut::next_base_segment()
{
    base_segment_i_ ++;
    if(base_segment_i_==base_segment_N_)base_segment_i_=0;
    show_base_segment();
}

void RobustCut::last_base_segment()
{
    if(base_segment_i_==0)base_segment_i_ = base_segment_N_ - 1;
    else base_segment_i_ --;
    show_base_segment();
}

void RobustCut::exit_base_segment()
{
    emit end();
}

void RobustCut::show_base_segment()
{
    std::vector<arma::umat>::iterator segiter;
    std::vector<arma::uvec>::iterator oiter;
    segiter = base_segment_list_.begin();
    oiter = labels_.begin();
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
//        std::cerr<<"segiter->n_cols:"<<segiter->n_cols<<std::endl;
        assert(base_segment_i_<segiter->n_cols);
        arma::uvec label = segiter->col(base_segment_i_);
        label += 1;
        mesh_ptr->graph_.sv2pix(label,*oiter);
        mesh_ptr->custom_color_.fromlabel(*oiter);
        ++oiter;
        if(oiter==labels_.end())break;
        ++segiter;
        if(segiter==base_segment_list_.end())break;
    }
    QString msg;
    emit message(msg.sprintf("Current:%u/%u,Next:->,Last:<-,Esc:Exit",base_segment_i_+1,base_segment_N_),0);
}

void RobustCut::consensus_segment()
{
    std::vector<arma::uvec>::iterator oiter;
    oiter = labels_.begin();
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    std::vector<arma::umat>::iterator segiter = base_segment_list_.begin();
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        arma::uvec label;
        cuts_.base_segments() = *segiter;
        cuts_.solve_consensus_segment(mesh_ptr,label);
        arma::uvec index = arma::find(*oiter==0);
        mesh_ptr->graph_.sv2pix(label,*oiter);
        if( !index.empty() && index.size() < oiter->size() )(*oiter)(index).fill(0);
        mesh_ptr->custom_color_.fromlabel(*oiter);
        ++oiter;
        if(oiter==labels_.end())break;
        ++segiter;
        if(segiter==base_segment_list_.end())break;
    }
    emit end();
}

