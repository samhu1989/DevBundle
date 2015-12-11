#ifndef QGLPOINTSELECTION_HPP
#define QGLPOINTSELECTION_HPP
#include "qglpointselection.h"
#include "common.h"
template<typename Mesh>
inline void PointSelections::selectAll(Mesh&m,arma::uvec&indices)
{
    for(iterator iter=begin();iter!=end();++iter)
    {
        arma::uvec result;
        selectAt(iter,m,result);
        indices = arma::join_cols(indices,result);
    }
    clear();
}

template<typename Mesh>
inline void PointSelections::selectAt(size_t index,Mesh&m,arma::uvec&indices)
{
    PointSelectionBase::Ptr ptr = (*this)[index];
    if(ptr&&0!=ptr.use_count())
    {
        switch(ptr->type())
        {
        case PointSelectionBase::RayPoint:
            {
                RayPointSelection::Ptr p = std::dynamic_pointer_cast<RayPointSelection>(ptr);
                p->select<Mesh>(m,indices);
            }
            break;
        default:
            std::cerr<<"Invalid Selection Type"<<std::endl;
        }
    }
}

template<typename Mesh>
inline void PointSelections::selectAt(PointSelections::iterator iter,Mesh&m,arma::uvec&indices)
{
    PointSelectionBase::Ptr ptr = *iter;
    if(ptr&&0!=ptr.use_count())
    {
        switch(ptr->type())
        {
        case PointSelectionBase::RayPoint:
            {
                RayPointSelection::Ptr p = std::dynamic_pointer_cast<RayPointSelection>(ptr);
                p->select<Mesh>(m,indices);
            }
            break;
        default:
            std::cerr<<"Invalid Selection Type"<<std::endl;
        }
    }
}


template<typename Mesh>
void RayPointSelection::select(Mesh& m,arma::uvec& indices )
{
    indices.reset();
    arma::fmat p((float*)m.points(),3,m.n_vertices(),true,true);
    arma::fvec dir = toward_ - from_;
    dir = arma::normalise(dir);
    p.each_col() -= from_;
    arma::frowvec distsa = arma::sum(arma::square(p));
    arma::frowvec distsb = arma::square( dir.t()*p );
    arma::frowvec dists = arma::sqrt( distsa - distsb );
    p*=-1.0;
    p = arma::normalise(p);
    arma::frowvec vars = dir.t()*p;
    arma::uvec within = arma::find( dists < 0.05 );
    if(within.is_empty())return;
    arma::uword select;
    vars(within).min(select);
    indices = arma::uvec(1);
    indices(0) = within(select);
}

#endif // QGLPOINTSELECTION_HPP

