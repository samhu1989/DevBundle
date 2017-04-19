#ifndef CUBE_H
#define CUBE_H
#include "common_global.h"
#include "MeshType.h"
#include <ext/hash_map>
#include <functional>
#include <QRgb>
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
    static Cube::Ptr newCube(void);
    static Cube::Ptr newCube(DefaultMesh&);
    static Cube::PtrLst newCubes(DefaultMesh& m, uint32_t N);
    static void colorByLabel(uint32_t* c, arma::uword size, const arma::uvec &label);
    static uint32_t colorFromLabel(uint32_t label);
    static void reset_color_set();
    void colorByLabel(uint32_t label);
    virtual void translate(
            const arma::fvec& t,
            Cube& result
            );
    virtual void transform(
            const arma::fmat& R,
            const arma::fvec& t,
            Cube& result
            );
    virtual void rotate(
            const arma::fmat& R,
            Cube& result
            );
    virtual void scale(
            const arma::fvec& s,
            Cube& result
            );
    virtual void scaleTo(
            const arma::fvec& s,
            Cube& result
            );
    virtual void scaleTo(
            const arma::fvec& s
            );
    virtual arma::vec get_dist2(const arma::fmat& v);
    virtual arma::vec get_dist2_for_plate(const arma::fmat& v,const arma::fvec& c, arma::uword zero_dim);
    virtual arma::vec dist(const arma::fmat&, const arma::fvec&, arma::uword zero_dim, arma::uword dim);
    virtual arma::vec get_dist2_box(const arma::fmat& v);
    arma::uvec outside(const arma::fmat& v);
    arma::uvec inside(const arma::fmat& v);

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
    virtual arma::fvec center_pos();
    virtual void updateScale();
    inline const arma::fvec& size()const{return size_;}
    arma::fmat corners_;
    arma::fmat weighted_corners_;
    arma::fvec obj_pos_;
    arma::uword label_;
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
    struct color_equal_to : public std::binary_function<uint32_t,uint32_t,bool>
    {
            bool operator()(const uint32_t& __x, const uint32_t& __y) const{
                float dr = ( qRed(__x) - qRed(__y) );
                float dg = ( qGreen(__x) - qGreen(__y) );
                float db = ( qBlue(__x) - qBlue(__y) );
                return std::sqrt( dr*dr + dg*dg + db*db ) < 50;
            }
    };
    static __gnu_cxx::hash_map<uint32_t,uint32_t,std::hash<uint32_t>,color_equal_to> color_label_;
    static __gnu_cxx::hash_map<uint32_t,uint32_t> label_color_;
};

}
#endif // CUBE_H
