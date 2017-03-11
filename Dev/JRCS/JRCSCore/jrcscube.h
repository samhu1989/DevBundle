#ifndef JRCSCUBE_H
#define JRCSCUBE_H
#include "jrcsbilateral.h"
namespace JRCS {
class JRCSCORESHARED_EXPORT Cube{
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
    arma::fmat corners_;
    arma::fmat weighted_corners_;
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
    arma::fvec obj_pos_;
    std::shared_ptr<arma::fmat> xv_;
    std::shared_ptr<arma::fmat> xn_;
    std::shared_ptr<arma::Mat<uint8_t>> xc_;

    std::shared_ptr<arma::fmat> param_mat_;
    std::shared_ptr<arma::fvec> param_vec_;
public:
    arma::fcube param_;
    const static uint32_t point_num_for_plate_ = 4;
    const static uint32_t plate_num_for_cube_ = 5;
    const static uint32_t point_num_for_cube_ = 20;
    static std::vector<arma::uvec> c4v_;
    static arma::fvec  scale_r_;
};
class JRCSCORESHARED_EXPORT JRCSCube:public JRCSBilateral
{
public:
    JRCSCube();
    virtual ~JRCSCube(){}
    virtual std::string name()const{return "JRCSCube";}
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//initialize the X
protected:
    virtual void compute(void);
    virtual void reset_obj_vn(
            float radius,
            arma::fvec& pos,
            arma::fmat& ov,
            arma::fmat& on
            );
protected:
    virtual void reset_alpha_cube();
    virtual void reset_prob_cube();
    virtual void prepare_cube();
    virtual void finish_cube();
    //calculate alpha
    //update r t
    //voting
    virtual void step_1(int i);
    virtual void updateRTforObj(
            const int start,
            const int end,
            arma::frowvec& colsum,
            arma::fmat& R,
            arma::fvec& t,
            Cube::PtrLst cube_ptrlst
            );
    //extracting new planes
    //updating var and p
    virtual void step_2(void);
    virtual bool isEnd_cube(void);
protected:
    Cube::PtrLst cube_ptrlst_;
    std::vector<Cube::PtrLst> cube_t_ptrlst_;
};
}
#endif // JRCSCUBE_H
