#ifndef MESHTYPE
#define MESHTYPE
#include <memory>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/StripifierT.hh>
#include "MeshColor.h"
struct Traits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<Traits>  DefaultMesh;

template <typename M>
class MeshBundle
{
public:
    using Ptr = std::shared_ptr<MeshBundle<M>>;
    MeshBundle():
        custom_color_(mesh_),
        strips_(mesh_)
    {}
    M                           mesh_;
    MeshColor<M>        custom_color_;
    OpenMesh::StripifierT<M>  strips_;
};

#endif // MESHTYPE

