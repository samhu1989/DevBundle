#ifndef DENSECRF3D_H
#define DENSECRF3D_H
#include "densecrf.h"

class SEGMENTATIONCORESHARED_EXPORT DenseCRF3D: public DenseCRF
{
public:
    DenseCRF3D(
            const arma::fmat& xyz,
            const arma::fmat& nxyz,
            const arma::fmat& rgb,
            const int L
            );
    virtual ~DenseCRF3D(){}
    void addPairwiseGaussian(
            arma::fvec sxyz,
            LabelCompatibility * function=NULL,
            KernelType kernel_type=DIAG_KERNEL,
            NormalizationType normalization_type=NORMALIZE_SYMMETRIC
            );
    void addPairwiseBilateral(
            arma::fvec sxyz,
            arma::fvec snxyz,
            arma::fvec srgb,
            LabelCompatibility * function=NULL,
            KernelType kernel_type=DIAG_KERNEL,
            NormalizationType normalization_type=NORMALIZE_SYMMETRIC
            );
    using DenseCRF::setUnaryEnergy;
protected:
    arma::fmat xyz_;
    arma::fmat nxyz_;
    arma::fmat rgb_;
};

#endif // DENSECRF3D_H
