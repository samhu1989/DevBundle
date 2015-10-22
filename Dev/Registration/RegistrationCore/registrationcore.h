#ifndef REGISTRATIONCORE_H
#define REGISTRATIONCORE_H
#include "registrationcore_global.h"
#include <OpenMesh/Core/Mesh/TriMeshT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <armadillo>
#include <iostream>
#include <memory>
namespace Registration
{
    namespace MeshPort{
        struct DefaultTraits : public OpenMesh::DefaultTraits
        {
            HalfedgeAttributes(OpenMesh::Attributes::PrevHalfedge);
        };
        typedef OpenMesh::TriMesh_ArrayKernelT<DefaultTraits>  DefaultMesh;
        template<typename MeshType>
        void meshToMat(MeshType&mesh,arma::fmat&X);
        template<typename MeshType>
        void matToMesh(arma::fmat&X,MeshType&mesh);
    }

}
#endif // REGISTRATIONCORE_H
