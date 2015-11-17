#ifndef OBJECTMODEL_H
#define OBJECTMODEL_H
#include <armadillo>
#include "common.h"
struct ObjModel
{
    typedef std::shared_ptr<ObjModel> Ptr;
    MeshBundle<DefaultMesh>::Ptr GeoM_;
    arma::fvec GeoVar_;
    arma::gmm_diag ColorM_;
};

#endif // OBJECTMODEL_H
