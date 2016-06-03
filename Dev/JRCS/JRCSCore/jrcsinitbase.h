#ifndef JRCSINITBASE_H
#define JRCSINITBASE_H
#include <memory>
#include <armadillo>
#include <common.h>
namespace JRCS
{
class JRCSInitBase
{
public:
    typedef std::shared_ptr<arma::fmat> MatPtr;
    typedef std::vector<MatPtr> MatPtrLst;
    typedef std::shared_ptr<arma::Mat<uint8_t>> CMatPtr;
    typedef std::vector<CMatPtr> CMatPtrLst;
    typedef std::shared_ptr<arma::Col<uint32_t>> LCMatPtr;
    typedef std::vector<LCMatPtr> LCMatPtrLst;
    typedef std::shared_ptr<arma::uvec> LMatPtr;
    typedef std::vector<LMatPtr> LMatPtrLst;
    JRCSInitBase();
    virtual bool configure(Config::Ptr config);
    virtual bool init_with_label(
            const int k,
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vlc,
            const LMatPtrLst& vl,
            int verbose
            );
    virtual void getAlpha(MatPtrLst&);
    virtual void getObjProb(arma::fvec&);

protected:
    bool check_centers();
    void extract_patch_features();
    void pca();
    void learn();
    void assign();
    void assign(const arma::mat& f,arma::fmat& p);
    void generate_alpha();
    void generate_prob();
private:
    //inputs:
    int k_;
    MatPtrLst vv_;
    MatPtrLst vn_;
    CMatPtrLst vc_;
    LCMatPtrLst vlc_;
    LMatPtrLst vl_;
    //output
    MatPtrLst alpha_;
    arma::fvec prob_;


    //clustering
    std::vector<arma::urowvec> input_patch_label_value_;
    std::vector<arma::mat> patch_features_;
    std::vector<arma::uvec> patch_sizes_;
    std::vector<arma::fmat> patch_prob_;
    arma::mat feature_base_;
    arma::mat feature_centers_;
    arma::gmm_diag gmm_;
    int verbose_;
    Config::Ptr config_;
};
}
#endif // JRCSINITBASE_H
