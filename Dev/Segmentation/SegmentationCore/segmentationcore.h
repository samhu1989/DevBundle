#ifndef SEGMENTATIONCORE_H
#define SEGMENTATIONCORE_H
#include "supervoxelclustering.h"
#include <armadillo>
#include <stdint.h>
#include "segmentationcore_global.h"
#include "MeshType.h"
namespace Segmentation
{
template class SEGMENTATIONCORESHARED_EXPORT  SuperVoxel<DefaultMesh>;
template class SEGMENTATIONCORESHARED_EXPORT  DefaultVoxelDistFunctor<DefaultMesh>;
template class SEGMENTATIONCORESHARED_EXPORT  SuperVoxelClustering<DefaultMesh>;
typedef SuperVoxelClustering<DefaultMesh> DefaultSuperVoxelClustering;
}

#endif // SEGMENTATIONCORE_H
