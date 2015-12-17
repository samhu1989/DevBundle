#ifndef FEATURECORE_H
#define FEATURECORE_H
#include "MeshType.h"
#include "featurecore_global.h"
#include <armadillo>
#include "pointnormal.h"
#include "colorhistogram.h"
#include "normalhistogram.h"
template class FEATURECORESHARED_EXPORT Feature::ColorHistogramLab<DefaultMesh>;
template class FEATURECORESHARED_EXPORT Feature::ColorHistogramRGB<DefaultMesh>;
template class FEATURECORESHARED_EXPORT Feature::NormalHistogram<DefaultMesh>;
#endif // FEATURECORE_H
