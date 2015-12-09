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
    void init(arma::fmat&X);
    void updateFullModel(MeshBundle<DefaultMesh>::Ptr);
    void updateColor(MeshBundle<DefaultMesh>::Ptr);
    void finishColor();
    void updateWeight(MeshBundle<DefaultMesh>::Ptr);
    void finishWeight();
    void computeLayout();
    void computeFullLayout();
    bool transform(DefaultMesh&,uint32_t);
    bool transform(DefaultMesh&,arma::fmat&,arma::fvec&);

    bool save(const std::string& path);
    bool load(const std::string& path);

    bool fullLayout(arma::fmat&,int32_t);
    bool fullModel(DefaultMesh&,int32_t);

    MeshBundle<DefaultMesh>::Ptr GeoM_;
    MeshBundle<DefaultMesh>::Ptr GeoLayout_;

    std::vector<T::Ptr> GeoT_;
    arma::fvec DistP_;
    arma::fvec ColorP_;
    arma::fvec NormP_;
private:
    DefaultMesh FullM_;
    DefaultMesh FullLayout_;
    arma::mat accu_color_;
    arma::uword accu_count_;
};

#endif // OBJECTMODEL_H
