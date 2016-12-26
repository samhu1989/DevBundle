#ifndef JRCSHOUGH_H
#define JRCSHOUGH_H
#include "jrcsbilateral.h"
#include <ext/hash_map>"
namespace JRCS{
struct  HoughSpace{
    void vote(
            const arma::fvec& v,
            const arma::fvec& n,
            const arma::Col<uint8_t>& c
            );
    void get_prob(
            const arma::fvec& v,
            const arma::fvec& n,
            const arma::Col<uint8_t>& c,
            arma::rowvec& prob
            );
    void extract(
            arma::fmat& v,
            arma::fmat& n,
            arma::Mat<uint8_t>& c
            );
private:
    arma::fmat R_;
    arma::fvec t_;
    //voting space
    arma::fcube param_z_u_;
    arma::fcube param_x_f_;
    arma::fcube param_x_b_;
    arma::fcube param_y_f_;
    arma::fcube param_y_b_;
    //color space
    arma::fvec  rgb_z_u_;
    arma::fvec  rgb_x_f_;
    arma::fvec  rgb_x_b_;
    arma::fvec  rgb_y_f_;
    arma::fvec  rgb_y_b_;
    //parameter
    arma::fvec  res_z_u_;
    arma::fvec  res_x_f_;
    arma::fvec  res_x_b_;
    arma::fvec  res_y_f_;
    arma::fvec  res_y_b_;
};
class JRCSCORESHARED_EXPORT JRCSHough:public JRCSBilateral
{
public:
    JRCSHough();
protected:
    //calculate alpha
    //update r t
    //voting
    void step_1(int i);
    //extracting new planes
    //updating var and p
    void step_2(void);
private:
    std::vector<std::shared_ptr<HoughSpace>> hough_ptrlst_;
};
}
#endif // JRCSHOUGH_H
