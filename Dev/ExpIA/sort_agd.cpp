#include "sort_agd.h"
#include "agd.h"
Sort_AGD::Sort_AGD(
        MeshBundle<DefaultMesh>::PtrList& inputs,
        QObject *parent
        ):QObject(parent),inputs_(inputs)
{
    ;
}

bool Sort_AGD::configure(Config::Ptr)
{
    if(inputs_.empty())
    {
        return false;
    }
    InputList::iterator iter;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>& m = **iter;
        if(m.graph_.empty())return false;
    }
    return true;
}

void Sort_AGD::process(void)
{
    InputList::iterator iter;
    Feature::AGD<DefaultMesh> agd;
    QString msg;
    uint32_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>& m = **iter;
        arma::vec agd_vec;
        emit message(msg.sprintf("Processing %u/%u",index,inputs_.size()),-1);
        agd.extract(m.graph_,agd_vec);
        sort(agd_vec,m);
        ++index;
        emit message(msg.sprintf("Processing %u/%u",index,inputs_.size()),-1);
    }
    emit message(tr("Sort_AGD::process() is finishing"),-1);
    emit finished();
}

void Sort_AGD::sort(const arma::vec& agd,MeshBundle<DefaultMesh>& m)
{
    std::cerr<<"sorting"<<std::endl;
//    std::cerr<<"agd.min():"<<agd.min()<<std::endl;
//    std::cerr<<"agd.max():"<<agd.max()<<std::endl;
//    std::cerr<<"agd:"<<std::endl;
//    std::cerr<<agd<<std::endl;
//    std::cerr<<"graph_.voxel_centers:"<<std::endl;
//    std::cerr<<m.graph_.voxel_centers<<std::endl;
//    std::cerr<<"graph_.voxel_neighbors:"<<std::endl;
//    std::cerr<<m.graph_.voxel_neighbors<<std::endl;
    arma::uvec descend_index = arma::sort_index(agd,"descend");
    arma::uvec descend_label(m.graph_.voxel_label.size(),arma::fill::zeros);
    arma::uvec::iterator iter_index;
    arma::uword current_process = 0;
    for( iter_index = descend_index.begin() ; iter_index != descend_index.end() ; ++iter_index )
    {
        arma::uvec current_index = arma::find( m.graph_.voxel_label == ( *iter_index + 1 ));
        descend_label.subvec(current_process,current_process+current_index.size()-1) = current_index;
        current_process += current_index.size();
    }
//    std::cerr<<"descend label:"<<std::endl;
//    std::cerr<<descend_label<<std::endl;
    arma::fmat v((float*)m.mesh_.points(),3,m.mesh_.n_vertices(),false,true);
//    std::cerr<<"v before:"<<std::endl;
//    std::cerr<<v<<std::endl;
    v = v.cols(descend_label);
//    std::cerr<<"v after:"<<std::endl;
//    std::cerr<<v<<std::endl;
    if(m.mesh_.has_vertex_normals())
    {
        arma::fmat n((float*)m.mesh_.vertex_normals(),3,m.mesh_.n_vertices(),false,true);
        n = n.cols(descend_label);
    }
    if(m.mesh_.has_vertex_colors())
    {
        arma::Mat<uint8_t> c((uint8_t*)m.mesh_.vertex_colors(),3,m.mesh_.n_vertices(),false,true);
        c = c.cols(descend_label);
    }
}
