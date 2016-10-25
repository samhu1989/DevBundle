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
    bool configure(Config::Ptr);
    void extract(const arma::mat& f,const arma::uvec& l,arma::mat& h);//for a frame
    arma::vec extract(const arma::mat& f);//for a patch
    void learn(const MatPtrLst& f, const LabelLst& l, MatPtrLst &h);
protected:
    void learn_code_book(const MatPtrLst& f);
private:
    arma::uword codebook_size_;
    arma::gmm_diag gmm_;
    arma::vec idf_;
};
}
#endif // BOF_H
