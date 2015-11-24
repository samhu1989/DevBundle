#ifndef VOXELGRAPH_H
#define VOXELGRAPH_H
#include "common.h"
class VoxelGraph
{
public:
    typedef std::shared_ptr<VoxelGraph> Ptr;
    typedef std::vector<Ptr> PtrList;
    VoxelGraph();
};

#endif // VOXELGRAPH_H
