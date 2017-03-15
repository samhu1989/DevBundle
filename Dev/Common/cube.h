#ifndef CUBE_H
#define CUBE_H
#include "common_global.h"
#include "MeshType.h"
namespace Common {
class COMMONSHARED_EXPORT Cube{
public:
    typedef std::shared_ptr<Cube> Ptr;
    typedef std::vector<Ptr> PtrLst;
    Cube();
    Cube(
        const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::fvec& pos
    );
    static Cube::Ptr newCube(DefaultMesh&);
    static Cube::PtrLst newCubes(DefaultMesh& m, uint32_t N);
    virtual void translate(
            const arma::fvec& t,
            Cube& result
            );
    virtual void transform(
            const arma::fmat& R,
            const arma::fvec& t,
            Cube& result
            );
    virtual void scale(
            const arma::fvec& s,
            Cube& result
            );
    virtual void scaleTo(
            const arma::fvec& s
            );
    virtual arma::vec get_dist2(const arma::fmat& v);
    virtual arma::vec get_dist2_for_plate(const arma::fmat& v,const arma::fvec& c, arma::uword zero_dim);
    virtual arma::vec dist(const arma::fmat&, const arma::fvec&, arma::uword zero_dim, arma::uword dim);
    virtual void get_weighted_corners(const arma::fmat& v, const arma::vec &alpha);
    virtual void get_weighted_color(const arma::fmat& v,const arma::Mat<uint8_t>& c );
    virtual void accumulate(
            const arma::fmat& v,
            const arma::fmat& n,
            const arma::Mat<uint8_t>& c,
            const arma::vec alpha
            );
    virtual void start_accumulate(const int r,const int c,const int s,const int num);
    virtual void accumulate(const Cube&,const int i);
    virtual void fit(void);
    virtual arma::fvec bottom_pos();
    virtual void updateScale();
    inline const arma::fvec& size()const{return size_;}
    arma::fmat corners_;
    arma::fmat weighted_corners_;
    arma::fvec obj_pos_;
protected:
    //use v to update plate centroids
    void median();
    void mean();
    void updateV2Centroids(void);
    void updateV2Corners(void);
    void updateCorners2V(void);
    void updateZeroDim(void);
    void updateCorners2Size(void);
private:
    arma::fvec size_;
    arma::fmat R_;
    arma::fvec t_;
    arma::fmat plate_centroids_;
    arma::uvec plate_zero_dim_;

    std::shared_ptr<arma::fmat> xv_;
    std::shared_ptr<arma::fmat> xn_;
    std::shared_ptr<arma::Mat<uint8_t>> xc_;

    std::shared_ptr<arma::fmat> param_mat_;
    std::shared_ptr<arma::fvec> param_vec_;
public:
    arma::fcube param_;
    const static uint32_t point_num_for_plate_;
    const static uint32_t plate_num_for_cube_;
    const static uint32_t point_num_for_cube_;
    static std::vector<arma::uvec> c4v_;
    static arma::fvec  scale_r_;
private:
    MeshBundle<DefaultMesh>::Ptr mesh_;
};

}
#endif // CUBE_H
