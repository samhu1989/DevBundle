#ifndef MESHTYPE
#define MESHTYPE
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
struct Traits : public OpenMesh::DefaultTraits
{
  HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
};

typedef OpenMesh::TriMesh_ArrayKernelT<Traits>  DefaultMesh;
#endif // MESHTYPE

