#include "gdcthread.h"
#include "meshpairviewerwidget.h"
#include "gdcoord.h"
#include <cassert>
GDCThread::GDCThread(
        std::vector<QWidget*> &inputs,
        QTimer &gl_timer,
        QObject* parent
        ):inputs_(inputs)
{
    ;
}

bool GDCThread::configure(Config::Ptr config)
{
    config_ = config;
    if(inputs_.empty())return false;
    int pick_num = -1;
    for(std::vector<QWidget*>::iterator iter = inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = qobject_cast<MeshPairViewerWidget*>(*iter);
        if( !w )return false;
        if( w->first_ptr()->graph_.empty() ){
            std::cerr<<"! Build Supervoxel First !"<<std::endl;
            return false;
        }
        if( w->first_selected().empty() )
        {
            std::cerr<<"! Please Select Points for Each Frame First !"<<std::endl;
            return false;
        }else if( pick_num < 0)
        {
            pick_num = w->first_selected().size();
        }else if( w->first_selected().size() != pick_num ){
            std::cerr<<"! Please Select Same Number of Points for Each Frame!"<<std::endl;
            return false;
        }
    }
    return true;
}

void GDCThread::process()
{
    for(std::vector<QWidget*>::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        MeshBundle<DefaultMesh>::Ptr m_ptr = w->first_ptr();
        arma::uvec selected_vox;
        m_ptr->graph_.getSvIndex(arma::uvec(w->first_selected()),selected_vox);
        arma::fmat vox_feature;
        Feature::GDCoord<DefaultMesh> f;
        f.extract(m_ptr->graph_,selected_vox,vox_feature);
        m_ptr->p_feature_ = arma::fmat(vox_feature.n_rows,m_ptr->mesh_.n_vertices());
        for(int r = 0 ; r < vox_feature.n_rows ; ++r )
        {
            arma::vec vf = arma::conv_to<arma::vec>::from(vox_feature.row(r).t());
            arma::vec pf;
            m_ptr->graph_.getPixFunc(vf,pf);
            m_ptr->p_feature_.row(r) = arma::conv_to<arma::frowvec>::from(pf.t());
        }
        feature_to_color((uint32_t*)m_ptr->custom_color_.vertex_colors(),m_ptr->mesh_.n_vertices(),m_ptr->p_feature_);
    }
    emit end();
}

void GDCThread::feature_to_color(uint32_t* ptr,arma::uword size,const arma::fmat& f)
{
    assert( f.n_cols == size );
    float max = f.max();
    float min = f.min();
//    std::cerr<<"f("<<min<<","<<max<<")"<<std::endl;
    #pragma omp parallel for
    for(int i = 0 ; i < f.n_cols ; ++i )
    {
        float r = 1.0 , g = 1.0 , b = 1.0 ;
        if(f.n_rows>0) r = ( f(0,i) - min ) / ( max - min );
        if(f.n_rows>1) g = ( f(1,i) - min ) / ( max - min );
        if(f.n_rows>2) b = ( f(2,i) - min ) / ( max - min );
        ptr[i] = qRgb( 255.0*r, 255.0*g, 255.0*b );
    }
}
