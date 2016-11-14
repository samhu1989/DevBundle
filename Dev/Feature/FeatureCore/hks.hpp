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
bool HKS<Mesh>::configure(Config::Ptr config)
{
    return true;
}
template<typename Mesh>
void HKS<Mesh>::extract(const typename MeshBundle<Mesh>::Ptr m,arma::mat& f)
{
    constructL(m);
    decomposeL();
    computeHKS(f);
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
        if( fiter == f.end() )break;
    }
}

template<typename Mesh>
void HKS<Mesh>::extract(const typename MeshBundle<Mesh>::Ptr m,const arma::uvec& l,arma::mat& f)
{
    f = arma::mat(6,m->graph_.size(),arma::fill::zeros);
    arma::uword label = 1;
    arma::uword max = arma::max(l);
    arma::uvec index;
    while( label <= max )
    {
        index = arma::find( label == l );
        if( index.size() < 5 ){
            ++label;
            if(label>max)break;
            else continue;
        }
//        std::cerr<<"constructL(m,index);"<<std::endl;
        constructL(m,index);
//        std::cerr<<"decomposeL();"<<std::endl;
        decomposeL();
        arma::mat patch_f;
//        std::cerr<<"computeHKS(patch_f);"<<std::endl;
        computeHKS(patch_f);
//        std::cerr<<"f.cols(index) = patch_f;"<<std::endl;
//        std::cerr<<"n:"<<index.size()<<std::endl;
//        std::cerr<<"n:"<<patch_f.n_cols<<std::endl;
        arma::uvec vox_index;
        m->graph_.getSvIndex(index,vox_index);
//        std::cerr<<"++"<<std::endl;
        assert(vox_index.size()==patch_f.n_cols);
        vox_index -= 1; // the voxel label is start from one
        f.cols(vox_index) = patch_f;
        std::cerr<<label<<"/"<<max<<std::endl;
        ++label;
//        std::cerr<<"__"<<std::endl;
    }
}

template<typename Mesh>
void HKS<Mesh>::extract(const typename MeshBundle<Mesh>::PtrList m,const LabelLst& l,MatPtrLst& f)
{
    f.resize(m.size());
    MatPtrLst::iterator fiter=f.begin();
    LabelLst::const_iterator liter=l.cbegin();
    for(typename MeshBundle<Mesh>::PtrList::const_iterator iter = m.begin() ; iter != m.end() ; ++iter )
    {
        (*fiter).reset(new arma::mat());
        extract(*iter,*liter,**fiter);
        ++ fiter;
        if( fiter == f.end() )break;
        ++ liter;
        if( liter == l.cend() )break;
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
void HKS<Mesh>::constructL(const typename MeshBundle<Mesh>::Ptr m,const arma::uvec& index)
{
    Mesh sub_mesh;
    typename VoxelGraph<Mesh>::Ptr graph_ptr = VoxelGraph<Mesh>::getSubGraphPtr(m->graph_,index,sub_mesh);
    const VoxelGraph<Mesh>& graph = *graph_ptr;
    size_t N = graph.size();
    W_ = arma::speye(N,N);
    for(arma::Mat<uint16_t>::const_iterator niter=graph.voxel_neighbors.begin();niter!=graph.voxel_neighbors.end();   )
    {
        uint16_t wi = *niter;
        ++niter;
        uint16_t wj = *niter;
        ++niter;
        assert(wi<graph.voxel_centers.n_cols);
        assert(wj<graph.voxel_centers.n_cols);
//        std::cerr<<"r:"<<graph.voxel_centers.n_rows<<std::endl;
//        std::cerr<<"c:"<<graph.voxel_centers.n_cols<<std::endl;
        double affinity = vecAffinity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_centers.col(wj),
                    d_scale_
                    );
//        std::cerr<<"ccc"<<std::endl;
        affinity = std::exp(affinity);
        affinity += convexity<arma::fvec>(
                    graph.voxel_centers.col(wi),
                    graph.voxel_normals.col(wi),
                    graph.voxel_centers.col(wj),
                    graph.voxel_normals.col(wj),
                    convex_scale_
                    );
//        std::cerr<<"ddd"<<std::endl;
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
}
template<typename Mesh>
void HKS<Mesh>::computeHKS(arma::mat&f)
{
    arma::mat sq_phi = arma::square(eig_vec_).t();
    arma::vec alpha_tao = arma::exp2(arma::linspace<arma::vec>(1.0,25.0,(25.0-1.0)*16.0+1.0));
    arma::mat HKS(alpha_tao.n_rows,W_.n_cols);
    #pragma omp for
    for(arma::uword i=0; i < HKS.n_rows ; ++i )
    {
        arma::vec lambda_a_tao = -lambda_*alpha_tao(i);
        HKS.row(i) = arma::sum(sq_phi.each_col()%arma::exp(lambda_a_tao));
    }
    //Logarithmic
    HKS = arma::log(HKS);
    //Derivative
    arma::mat Dif_HKS(HKS.n_rows-1,HKS.n_cols);
    #pragma omp for
    for(arma::uword i=0; i < Dif_HKS.n_rows ; ++i )
    {
        Dif_HKS.row(i) = HKS.row(i+1) - HKS.row(i);
    }
    //Fourier Transform
    arma::mat SI_HKS = arma::abs(arma::fft(Dif_HKS));
    f = SI_HKS.rows(0,5);
}
}
