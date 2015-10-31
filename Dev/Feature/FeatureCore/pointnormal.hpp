#include "pointnormal.h"
#include "nanoflann.hpp"
using namespace nanoflann;
namespace Feature{
template<typename M>
void computePointNormal(M& mesh,float radius,int k)
{
    MeshKDTreeInterface<M> points(mesh);
    KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<float,MeshKDTreeInterface<M>>,
            MeshKDTreeInterface<M>,
            3>
            kdtree(3,points,KDTreeSingleIndexAdaptorParams(9));
    kdtree.buildIndex();
    nanoflann::SearchParams param;

}
}

