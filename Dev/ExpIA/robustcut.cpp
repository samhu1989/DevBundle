#include "robustcut.h"
#include <QMessageBox>
RobustCut::RobustCut(
        MeshBundle<DefaultMesh>::PtrList &inputs,
        std::vector<arma::uvec> &labels,
        QObject *parent
        ):inputs_(inputs),labels_(labels),QObject(parent)
{
    setObjectName("RobustCut");
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
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {

        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        arma::uvec label;
        cuts_.cutGraph(mesh_ptr,label);
        base_segment_list_.emplace_back(cuts_.base_segments());
        arma::uword n = base_segment_list_.back().n_cols;
        if(n<base_segment_N_)base_segment_N_=n;
        ++oiter;
        if(oiter==labels_.end())break;
    }
    emit message(tr("Next:->,Last:<-,Esc:Exit"),0);
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
        assert(base_segment_i_<segiter->n_cols);
        arma::uvec label = segiter->col(base_segment_i_);
        mesh_ptr->graph_.sv2pix(label,*oiter);
        mesh_ptr->custom_color_.fromlabel(*oiter);
        ++oiter;
        if(oiter==labels_.end())break;
        ++segiter;
        if(segiter==base_segment_list_.end())break;
    }
    emit message(tr("Next:->,Last:<-,Esc:Exit"),0);
}

void RobustCut::consensus_segment()
{
    std::vector<arma::uvec>::iterator oiter;
    oiter = labels_.begin();
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {

        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        arma::uvec label;
        cuts_.cutGraph(mesh_ptr,label);
        mesh_ptr->graph_.sv2pix(label,*oiter);
        mesh_ptr->custom_color_.fromlabel(*oiter);
        ++oiter;
        if(oiter==labels_.end())break;
    }
    emit end();
}

