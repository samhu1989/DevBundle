#ifndef BOF_H
#define BOF_H
#include "common.h"
#include "featurecore.h"
/***********************
 * Bag of Feature
 */
namespace Feature {
class FEATURECORESHARED_EXPORT BOF
{
public:
    typedef std::shared_ptr<arma::mat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    typedef std::vector<arma::uvec> LabelLst;
    typedef std::vector<arma::vec> VecLst;
    BOF();
    bool configure(Config::Ptr);
    inline void set_size(const arma::uword& size){codebook_size_=size;}
    void extract(const arma::mat& f,const arma::uvec& l,arma::mat& h);//for a frame
    void learn(const MatPtrLst& f, const LabelLst& l, MatPtrLst &h);
    inline const arma::mat gmm_mean(void)const{ return gmm_.means;}
    inline const VecLst& idf(void)const{return idf_;}
    inline const std::vector<arma::uvec>& assignment()const{return assignment_;}
protected:
    void learn_code_book(const MatPtrLst& f);
private:
    arma::uword codebook_size_;
    arma::gmm_diag gmm_;
    VecLst idf_;
    arma::vec g_idf_;
    std::vector<arma::uvec> assignment_;
    Config::Ptr config_;
};
}
#endif // BOF_H
