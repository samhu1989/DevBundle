#ifndef HIERARCHICALIZATION_H
#define HIERARCHICALIZATION_H
#include "common.h"
#include "segmentationcore_global.h"
#include <armadillo>
#include <nanoflann.hpp>
typedef struct
{
    arma::fvec::fixed<3> normal;
    std::vector<arma::uword> index;	// The neighborhood
    arma::fvec::fixed<3> eval;
    bool candidate;
}Neighbor;

typedef struct
{
    arma::fmat boxmat;
    arma::fvec center;
    float width;
    float height;
    float depth;
}BBox;

typedef struct
{
    float size_x;
    float size_y;
    float size_z;
    arma::uword id_;
    arma::uword father_id_;
    std::vector<arma::uword> child_plane_;
    std::vector<arma::uword> child_obj_;
    BBox bbox_;
    bool is_obj_;
}IdNode;

class SEGMENTATIONCORESHARED_EXPORT Hierarchicalization
{
public:
    typedef
    nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword>
    KDTree;
    Hierarchicalization();
    bool configure(Config::Ptr config);
    void compute(DefaultMesh&);
    void getObjectLabel(arma::uvec&);
    void getBBoxLabel(arma::uvec&);
    void getObjectBox(DefaultMesh&);
    void getPlaneLabel(arma::uvec&);
protected:
    void reset(const DefaultMesh&);
    void build(DefaultMesh&);
    void calsize(const arma::fmat& cloud,const arma::uword);
    void calneighbor(DefaultMesh&);
    float calarea(float size_x, float size_y, float size_z);
    BBox calbbox(arma::fmat&);
    float angle(const arma::fvec&,const arma::fvec&);
    void regiongrow(const arma::fmat& cloud,const arma::uword);
    arma::uvec withinNode(arma::uword);
    arma::uvec withinBox(const arma::fmat&);
private:
    uint32_t iden;
    uint32_t planenum;
    std::vector<IdNode> idtree_;
    std::vector<Neighbor> nei;	//output, every voxel's neighborhood
    arma::fmat cloud_;
    arma::ivec label_;
    bool  force_new_normal_;
    float neighbor_radius_;
    float anglethres_tight_;
    float anglethres_relax_;
    float point2plane_th_;
};


#endif // HIERARCHICALIZATION_H
