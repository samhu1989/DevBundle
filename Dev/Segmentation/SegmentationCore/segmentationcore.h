#ifndef SEGMENTATIONCORE_H
#define SEGMENTATIONCORE_H
#include "supervoxelclustering.h"
#include <armadillo>
#include <stdint.h>
#include "segmentationcore_global.h"
namespace Segmentation
{
template class SEGMENTATIONCORESHARED_EXPORT  SuperVoxel<DefaultMesh>;
template class SEGMENTATIONCORESHARED_EXPORT  DefaultVoxelDistFunctor<DefaultMesh>;
template class SEGMENTATIONCORESHARED_EXPORT  SuperVoxelClustering<DefaultMesh>;
typedef SuperVoxelClustering<DefaultMesh> SuperVox;
}

#endif // SEGMENTATIONCORE_H
