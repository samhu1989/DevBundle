#include "hks.h"
#include <armadillo>
namespace Feature{
template<typename Mesh>
HKS<Mesh>::HKS()
{
    d_scale_ = 9e-4;
    convex_scale_ = 0.01;
}
template<typename Mesh>
void HKS<Mesh>::extract(const typename MeshBundle<Mesh>::Ptr m,arma::mat& f)
{
    constructL(m);
    decomposeL();
    computeHKS();
}
template<typename Mesh>
void HKS<Mesh>::extract(const typename MeshBundle<Mesh>::PtrList m,MatPtrLst& f)
{
    f.resize(m.size());
    MatPtrLst::iterator fiter=f.begin();
    for(typename MeshBundle<Mesh>::PtrList::const_iterator iter = m.begin() ; iter != m.end() ; ++iter )
    {
        (*fiter).reset(new arma::mat());
        extract(*iter,**fiter);
        ++ fiter;
        if( fiter != f.end() )break;
    }
}
template<typename Mesh>
void HKS<Mesh>::constructL(const typename MeshBundle<Mesh>::Ptr m)
{
    const MeshBundle<Mesh>& mesh = *m;
    const VoxelGraph<Mesh>& graph = mesh.graph_;
    size_t N = graph.size();
    W_ = arma::speye(N,N);
    for(arma::Mat<uint16_t>::const_iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        double affinity = vecAffinity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_centers.col(wj),
                    d_scale_
                    );
        affinity = std::exp(affinity);
        affinity += convexity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_normals.col(wi),
                    graph.voxel_centers.col(wj),
                    graph.voxel_normals.col(wj),
                    convex_scale_
                    );
        W_(wi,wj) = 0.5*affinity;
        W_(wj,wi) = W_(wi,wj);
    }
    arma::vec D = arma::vectorise(arma::mat(arma::sum(W_)));
    arma::sp_mat Dmat = arma::speye<arma::sp_mat>(W_.n_rows,W_.n_cols);
    Dmat.diag() = D;
    W_ *= -1.0;
    W_ += Dmat;
}
template<typename Mesh>
void HKS<Mesh>::decomposeL()
{
    arma::uword k = 200;
    k = std::min(k,W_.n_cols-1);
    bool success = false;
    double stol = 50.0;
    double etol = std::numeric_limits<double>::epsilon();
    success = arma_custom::eigs_sym(lambda_,eig_vec_,W_,k,"sm",stol,etol);
    if(!success)std::cerr<<"Failed on decomposition, Please relax the tol"<<std::endl;//failed
//    arma::uvec index = arma::find( lambda_ <= 0.0 );
//    if(!index.empty()){
//        lambda_.shed_rows(index(0),index(index.size()-1));
//        eig_vec_.shed_cols(index(0),index(index.size()-1));
//    }
}
template<typename Mesh>
void HKS<Mesh>::computeHKS()
{
    ;
}
}
