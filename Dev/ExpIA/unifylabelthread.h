#ifndef UNIFYLABELCOLORSIZETHREAD_H
#define UNIFYLABELCOLORSIZETHREAD_H
#include <QThread>
#include "common.h"
#include <armadillo>
class UnifyLabelThread:public QThread
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator InputIterator;
    typedef std::vector<arma::uvec>::iterator OutputIterator;
    UnifyLabelThread(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& outputs,
            arma::mat& base,
            arma::mat& center
            ):inputs_(inputs),labels_(outputs),feature_base_(base),feature_centers_(center)
    {
        setObjectName("UnifyLabelThread");
    }
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
protected:
    void run();
    void extract_patch_features();
    void learn();
    void assign();
    void assign(const arma::mat& features,arma::urowvec& label_value);
    void alter_label(
            const arma::urowvec&,
            const arma::urowvec&,
            arma::uvec&
            );


private:
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<arma::uvec>& labels_;
    arma::mat& feature_base_;
    arma::mat& feature_centers_;
    std::vector<arma::mat> patch_features_;//each column is a feature vector for a patch
    std::vector<arma::urowvec> input_patch_label_value_;
    Config::Ptr config_;
    arma::gmm_diag gmm_;
};

#endif // UNIFYLABELCOLORSIZETHREAD_H
