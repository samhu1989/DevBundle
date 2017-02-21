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
    void translate(
            const arma::fvec& t,
            Plate& result
            );
    void get_local_translate(
            arma::fvec& t
            );
    void local_translate(
            const arma::fvec& t,
            Plate& result
            );
    void transform(
            const arma::fmat& R,
            const arma::fvec& t,
            Plate& result
            );
    void scale(
            const arma::fvec& s,
            Plate& result
            );
    arma::vec get_dist2(const arma::fmat& v);
    arma::vec dist(const arma::fmat& v,arma::uword dim);
    void get_weighted_centroid(const arma::fmat& v, const arma::vec &alpha);
    void accumulate(
            const arma::fmat& v,
            const arma::fmat& n,
            const arma::Mat<uint8_t>& c,
            const arma::vec alpha
            );
    void accumulate(const Plate&);
    void fit(void);
    void print(void);
    double area(void);
    arma::fvec size_;
    arma::fmat R_;
    arma::fvec t_;
    arma::fmat corners_;
    arma::fvec centroid_;
    arma::fvec origin_;
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
    void step_1(int i);
    //extracting new planes
    //updating var and p
    void step_2(void);
    virtual bool isEnd_primitive(void);
private:
    Plate::PtrLst plate_ptrlst_;
    std::vector<Plate::PtrLst> plate_t_ptrlst_;
    int plate_num_for_obj_;
    int point_num_for_plate_;
};
}
#endif // JRCSHOUGH_H
