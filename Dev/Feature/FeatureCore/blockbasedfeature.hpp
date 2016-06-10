#include "blockbasedfeature.h"
#include "common.h"
namespace Feature {
template<typename Mesh>
BlockBasedFeature<Mesh>::BlockBasedFeature()
{

}
template<typename Mesh>
void BlockBasedFeature<Mesh>::extract(const Mesh& mesh, arma::vec& feature)
{
    feature = arma::vec(12,arma::fill::zeros);
    size_t N = mesh.n_vertices();
    assert(N>2);
    arma::fmat points((float*)mesh.points(),3,N,false,true);
    arma::fmat box;
    get3DMBB(points,2,box);
    //length width height
    feature(0) = arma::norm( box.col(0) - box.col(1) );
    feature(1) = arma::norm( box.col(0) - box.col(3) );
    feature(2) = arma::norm( box.col(0) - box.col(4) );
    //six plane coefficient
    arma::fvec up;
    arma::uvec upIdx = {4,5,6,7};
    arma::fmat upV = box.cols(upIdx);
    fitPlane(upV,up);
    arma::fvec down;
    arma::uvec downIdx = {0,1,2,3};
    arma::fmat downV = box.cols(downIdx);
    fitPlane(downV,down);
    arma::fvec front;
    arma::uvec frontIdx = {0,1,5,4};
    arma::fmat frontV = box.cols(frontIdx);
    fitPlane(frontV,front);
    arma::fvec back;
    arma::uvec backIdx = {3,7,6,2};
    arma::fmat backV = box.cols(backIdx);
    fitPlane(backV,back);
    arma::fvec left;
    arma::uvec leftIdx = {1,2,6,5};
    arma::fmat leftV = box.cols(leftIdx);
    fitPlane(leftV,left);
    arma::fvec right;
    arma::uvec rightIdx = {0,4,7,3};
    arma::fmat rightV = box.cols(rightIdx);
    fitPlane(rightV,right);
    if( feature(0)<feature(1) )
    {
        feature.swap_rows(0,1);
        arma::fvec tmp = front;
        front = right;
        right = back;
        back = left;
        left = tmp;
    }
    //front back
    //left right
    //up down
    arma::mat dists(3,points.n_cols);
    for(size_t idx=0;idx<points.n_cols;++idx)
    {
        arma::fvec p(4);
        arma::uword minidx;
        p.head(3) = points.col( idx );
        p(3) = 1.0;
        dists(0,idx) = std::min(std::abs(arma::dot(p,front)),std::abs(arma::dot(p,back)));
        dists(1,idx) = std::min(std::abs(arma::dot(p,left)),std::abs(arma::dot(p,right)));
        dists(2,idx) = std::min(std::abs(arma::dot(p,up)),std::abs(arma::dot(p,down)));
        dists.col(idx).min(minidx);
        feature( 3 + minidx ) += 1.0;
        feature( 6 + minidx ) += dists( minidx , idx );
    }
    if(feature(3)>0)feature(6) /= feature(3);
    if(feature(4)>0)feature(7) /= feature(4);
    if(feature(5)>0)feature(8) /= feature(5);
    for(size_t idx=0;idx<points.n_cols;++idx)
    {
       arma::uword minidx;
       double mindist = dists.col(idx).min(minidx);
       feature( 9 + minidx ) += std::pow( mindist - feature( 6 + minidx ) , 2.0 );
    }
    if(feature(3)>1)
    {
        feature(9) /= ( feature(3) - 1.0 );
        feature(9) = std::sqrt(feature(9));
    }else {feature(9)=0;}
    if(feature(4)>1)
    {
        feature(10) /= ( feature(4) - 1.0 );
        feature(10) = std::sqrt(feature(10));
    }else {feature(10)=0;}
    if(feature(5)>1)
    {
        feature(11) /= ( feature(5) - 1.0 );
        feature(11) = std::sqrt(feature(11));
    }else {feature(11)=0;}
    feature(3) /= N;
    feature(4) /= N;
    feature(5) /= N;
    if(!feature.is_finite())
    {
        std::cerr<<"infinite block feature"<<std::endl;
        std::cerr<<feature<<std::endl;
    }
    assert(feature.is_finite());
}
}

