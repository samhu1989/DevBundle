#ifndef OBJECTMODEL_H
#define OBJECTMODEL_H
#include <armadillo>
#include "common.h"
struct ObjModel
{
    typedef std::shared_ptr<ObjModel> Ptr;
    typedef struct T{
      typedef std::shared_ptr<struct T> Ptr;
      float R[9];
      float t[3];
    }T;
    ObjModel();
    void update(MeshBundle<DefaultMesh>::Ptr ptr);
    void computeLayout();
    bool save(const std::string& path);
    bool load(const std::string& path);
    MeshBundle<DefaultMesh>::Ptr GeoM_;
    MeshBundle<DefaultMesh>::Ptr GeoLayout_;
    std::vector<T::Ptr> GeoT_;
    std::vector<float> GeoP_;
    std::shared_ptr<arma::gmm_diag> ColorM_;//not used yet
    std::vector<float> ColorP_;
};

#endif // OBJECTMODEL_H
