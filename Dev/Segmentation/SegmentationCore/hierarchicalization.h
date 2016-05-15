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
    void getObjectLabel(arma::uvec&lbl);
    void getPlaneLabel(arma::uvec&lbl);
    void getNonPlaneLabel(arma::uvec&lbl);
protected:
    void reset(const DefaultMesh&);
    void build(DefaultMesh&);
    void calsize(const arma::fmat& cloud,const arma::uword&);
    void calneighbor(DefaultMesh&);
    float calarea(float size_x, float size_y, float size_z);
    BBox calbbox(arma::fmat&);
    float angle(const arma::fvec&,const arma::fvec&);
    void regiongrow(const arma::fmat& cloud,const arma::uword&);
    void planegrow(const arma::uword&);
private:
    uint32_t iden;
    uint32_t planenum;
    std::vector<IdNode> idtree_;
    std::vector<Neighbor> nei;	//output, every voxel's neighborhood
    arma::ivec label_;
    float neighbor_radius_;
    float anglethres_;
};


#endif // HIERARCHICALIZATION_H
