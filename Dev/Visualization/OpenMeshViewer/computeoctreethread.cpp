#include "computeoctreethread.h"
void ComputeOctreeThread::run()
{
    DefaultOctree octree;
    DefaultOctreeContainer container(input_);
    octree.initialize(container,unibn::OctreeParams(float(0.5)));
    OctreeMeshInterface<DefaultOctree> oct_interface(octree,output_);
    oct_interface.getLeafMesh(true);
}

