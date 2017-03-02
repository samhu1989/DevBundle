#ifndef JRCSPRIMITIVE_H
#define JRCSPRIMITIVE_H
#include "jrcsbilateral.h"
#include <ext/hash_map>"
namespace JRCS{
struct  JRCSCORESHARED_EXPORT Plate{
    typedef std::shared_ptr<Plate> Ptr;
    typedef std::vector<Ptr> PtrLst;
//    /*   Y axis
//     *   ^
//     *   |
//     *   B
//     * C E A -> X axis
//     *   D
//     *
//     */
//    typedef enum{
//        A,B,C,D,E
//    }TYPE;
    Plate();
    Plate(const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c,
        const arma::fvec &pos
          );
    virtual void translate(
            const arma::fvec& t,
            Plate& result
            );
    virtual void get_local_translate(
            arma::fvec& t
            );
    virtual void local_translate(
            const arma::fvec& t,
            Plate& result
            );
    virtual void transform(
            const arma::fmat& R,
            const arma::fvec& t,
            Plate& result
            );
    virtual void scale(
            const arma::fvec& s,
            Plate& result
            );
    virtual arma::vec get_dist2(const arma::fmat& v);
    virtual arma::vec dist(const arma::fmat&,const arma::fvec&,arma::uword dim);
    virtual void get_weighted_centroid(const arma::fmat& v, const arma::vec &alpha);
    virtual void get_weighted_color(const arma::fmat& v,const arma::Mat<uint8_t>& c );
    virtual void accumulate(
            const arma::fmat& v,
            const arma::fmat& n,
            const arma::Mat<uint8_t>& c,
            const arma::vec alpha
            );
    virtual void start_accumulate(const int r,const int c,const int s,const int num);
    virtual void accumulate(const Plate&,const int i);
    virtual void median();
    virtual void mean();
    virtual void norm_mean();
    virtual void fit(void);
    virtual void print(void);
    virtual double area(void);
    arma::fvec size_;
    arma::fmat R_;
    arma::fvec t_;
    arma::fmat corners_;
    arma::fvec centroid_;
//    arma::fvec origin_;
    arma::fvec weighted_centroid_;
    arma::fvec obj_pos_;
    std::shared_ptr<arma::fmat> xv_;
    std::shared_ptr<arma::fmat> xn_;
    std::shared_ptr<arma::Mat<uint8_t>> xc_;
    //scale dim0
    //scale dim1
    //dt    dim2
    arma::fvec  scale_r_;
    arma::fvec  trans_r_;
    arma::fcube param_;
    std::shared_ptr<arma::fmat> param_mat_;
    std::shared_ptr<arma::fvec> param_vec_;
//    TYPE type_;
};

class JRCSCORESHARED_EXPORT JRCSPrimitive:public JRCSBilateral
{
public:
    JRCSPrimitive();
    virtual ~JRCSPrimitive(){}
    virtual std::string name()const{return "JRCSPrimitive";}
    virtual void initx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );//randomly initialize the X
protected:
    virtual void reset_obj_vn(
            float radius,
            arma::fvec& pos,
            arma::fmat& ov,
            arma::fmat& on
            );
    virtual void compute(void);
protected:
    virtual void reset_alpha_primitive();
    virtual void reset_prob_primitive();
    virtual void prepare_primitive();
    virtual void finish_primitive();
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
            Plate::PtrLst plate_ptrlst
            );
    //extracting new planes
    //updating var and p
    virtual void step_2(void);
    virtual bool isEnd_primitive(void);
protected:
    Plate::PtrLst plate_ptrlst_;
    std::vector<Plate::PtrLst> plate_t_ptrlst_;
    int plate_num_for_obj_;
    int point_num_for_plate_;
};
}
#endif // JRCSHOUGH_H
