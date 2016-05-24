#ifndef FEATURECORE_H
#define FEATURECORE_H
#include "MeshType.h"
#include "featurecore_global.h"
#include <armadillo>
#include "pointnormal.h"
#include "colorhistogram.h"
#include "normalhistogram.h"
#include "blockbasedfeature.h"
#include "blockbasedfeature.hpp"
template class FEATURECORESHARED_EXPORT Feature::ColorHistogramLab<DefaultMesh>;
template class FEATURECORESHARED_EXPORT Feature::ColorHistogramRGB<DefaultMesh>;
template class FEATURECORESHARED_EXPORT Feature::NormalHistogram<DefaultMesh>;
template class FEATURECORESHARED_EXPORT Feature::BlockBasedFeature<DefaultMesh>;
void FEATURECORESHARED_EXPORT extract_patch_feature(DefaultMesh&,arma::vec&,Config::Ptr);
#endif // FEATURECORE_H
