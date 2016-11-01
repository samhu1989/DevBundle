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
    BOF();
    inline void set_size(const arma::uword& size){codebook_size_=size;}
    void extract(const arma::mat& f,const arma::uvec& l,arma::mat& h);//for a frame
    arma::vec extract(const arma::mat& f);//for a patch
    void learn(const MatPtrLst& f, const LabelLst& l, MatPtrLst &h);
    const arma::mat gmm_mean(void)const{ return gmm_.means;}
protected:
    void learn_code_book(const MatPtrLst& f);
private:
    arma::uword codebook_size_;
    arma::gmm_diag gmm_;
    arma::vec idf_;
};
}
#endif // BOF_H
