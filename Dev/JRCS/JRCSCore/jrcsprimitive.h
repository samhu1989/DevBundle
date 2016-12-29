#ifndef JRCSPRIMITIVE_H
#define JRCSPRIMITIVE_H
#include "jrcsbilateral.h"
#include <ext/hash_map>"
namespace JRCS{
struct  Plate{
    typedef std::shared_ptr<Plate> Ptr;
    typedef std::vector<Ptr> PtrLst;
    Plate();
    Plate(
        const arma::fmat& v,
        const arma::fmat& n,
        const arma::Mat<uint8_t>& c
    );
    void transform(
            const arma::fmat& R,
            const arma::fvec& t,
            Plate&
            );
    arma::vec get_alpha(const arma::fmat& v);
    std::shared_ptr<arma::fmat> xv_;
    std::shared_ptr<arma::fmat> xn_;
    std::shared_ptr<arma::Mat<uint8_t>> xc_;
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
