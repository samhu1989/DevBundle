#ifndef MESHTYPE_H
#define MESHTYPE_H
#include <memory>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include "MeshColor.h"
struct Traits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};


typedef OpenMesh::PolyMesh_ArrayKernelT<Traits>  DefaultMesh;

template <typename M>
class MeshBundle
{
public:
    typedef typename std::shared_ptr<MeshBundle<M>> Ptr;
    MeshBundle():
        custom_color_(mesh_),
        strips_(mesh_)
    {}
    M                           mesh_;
    MeshColor<M>        custom_color_;
    OpenMesh::StripifierT<M>  strips_;
};

#endif // MESHTYPE

