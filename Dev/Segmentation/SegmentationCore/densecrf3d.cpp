#include "densecrf3d.h"

DenseCRF3D::DenseCRF3D(
        const arma::fmat& xyz,
        const arma::fmat& nxyz,
        const arma::fmat& rgb,
        const int L
        ):
    DenseCRF(xyz.n_cols,L),
    xyz_(xyz),
    nxyz_(nxyz),
    rgb_(rgb)
{
    ;
}

void DenseCRF3D::addPairwiseGaussian(
        arma::fvec sxyz,
        LabelCompatibility * function,
        KernelType kernel_type,
        NormalizationType normalization_type
        )
{
    arma::fmat feature(3,xyz_.n_cols,arma::fill::zeros);
    #pragma omp for
    for(int idx=0;idx<xyz_.n_cols;++idx)
    {
        feature(0,idx) = xyz_(0,idx) / sxyz(0);
        feature(1,idx) = xyz_(1,idx) / sxyz(1);
        feature(2,idx) = xyz_(2,idx) / sxyz(2);
    }
    addPairwiseEnergy( arma::conv_to<arma::mat>::from(feature), function, kernel_type, normalization_type );
}

void DenseCRF3D::addPairwiseBilateral(
        arma::fvec sxyz,
        arma::fvec snxyz,
        arma::fvec srgb,
        LabelCompatibility * function,
        KernelType kernel_type,
        NormalizationType normalization_type
        )
{
    arma::fmat feature(6,xyz_.n_cols,arma::fill::zeros);
    #pragma omp for
    for(int idx = 0 ; idx < xyz_.n_cols ; ++idx)
    {
        feature(0,idx) = xyz_(0,idx) / sxyz(0);
        feature(1,idx) = xyz_(1,idx) / sxyz(1);
        feature(2,idx) = xyz_(2,idx) / sxyz(2);
        feature(3,idx) = rgb_(0,idx) / srgb(0);
        feature(4,idx) = rgb_(1,idx) / srgb(1);
        feature(5,idx) = rgb_(2,idx) / srgb(2);
//        feature(6,idx) = nxyz_(0,idx) / snxyz(0);
//        feature(7,idx) = nxyz_(1,idx) / snxyz(1);
//        feature(8,idx) = nxyz_(2,idx) / snxyz(2);
    }
    addPairwiseEnergy( arma::conv_to<arma::mat>::from(feature), function, kernel_type, normalization_type );
}
