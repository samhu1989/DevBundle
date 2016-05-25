#ifndef JRCSINITBASE_H
#define JRCSINITBASE_H
#include <memory>
#include <armadillo>
#include <common.h>
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
    bool configure(Config::Ptr config);
    void init_alpha_with_label(
            MatPtrLst& alpha,
            const MatPtrLst& vv,
            const MatPtrLst& vn,
            const CMatPtrLst& vc,
            const LCMatPtrLst& vlc,
            const LMatPtrLst& vl,
            bool verbose
            );
protected:
    bool check_centers();
    void extract_patch_features();
    void pca();
    void learn();
    void generate_alpha();
private:
    //inputs:
    MatPtrLst vv_;
    MatPtrLst vn_;
    CMatPtrLst vc_;
    LCMatPtrLst vlc_;
    LMatPtrLst vl_;
    //
    std::vector<arma::urowvec> input_patch_label_value_;
    std::vector<arma::mat> patch_features_;
    arma::mat feature_base_;
    arma::mat feature_centers_;
    arma::gmm_diag gmm_;
    Config::Ptr config_;
};

#endif // JRCSINITBASE_H
