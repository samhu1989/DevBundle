#ifndef OCTREEGRID_H
#define OCTREEGRID_H
#include "filtercore.h"
namespace Filter {
template<typename M>
class OctreeGrid:public FilterBase<M>
{
public:
    typedef MeshOctreeContainer<M> OctreeMesh;
    typedef unibn::Octree<arma::fvec,MeshOctreeContainer<M>> Octree;
    typedef typename M::VertexHandle MeshVertexHandle;
    typedef typename M::Point MeshPoint;
    typedef typename M::Color MeshColor;
    OctreeGrid():
        FilterBase<M>(),
        seed_resolution_(0.0)
    {}
    inline void set_seed_resolution(double res){seed_resolution_=res;}
protected:
    virtual bool initCompute();
    virtual void applyFilter(M&);
private:
    double seed_resolution_;
};
}
#include "octreegrid.hpp"
template class FILTERCORESHARED_EXPORT Filter::OctreeGrid<DefaultMesh>;
#endif // OCTREEGRID_H
