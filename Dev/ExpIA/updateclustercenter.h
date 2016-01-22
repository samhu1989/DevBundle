#ifndef UPDATECLUSTERCENTER_H
#define UPDATECLUSTERCENTER_H
#include <QThread>
#include "common.h"
#include "objectmodel.h"
class UpdateClusterCenter : public QThread
{
    Q_OBJECT
public:
    typedef std::shared_ptr<DefaultMesh> MeshPtr;
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> InputList;
    typedef std::vector<arma::uvec> LabelList;
    typedef std::vector<ObjModel::Ptr> ObjList;
    UpdateClusterCenter(
            InputList& inputs,
            LabelList& labels,
            ObjList& objects,
            arma::mat& base,
            arma::mat& center
            ):inputs_(inputs),
              labels_(labels),
              objects_(objects),
              feature_base_(base),
              feature_centers_(center),
              raw_feature_dim(0)
    {
        setObjectName("UpdateClusterCenter");
    }
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
protected:
    void run();
    void evaluate_patches();
    double evaluate_patch(
            uint64_t oidx,
            uint64_t fidx,
            DefaultMesh& patch
            );
    double match_patch(
            DefaultMesh& om,
            DefaultMesh& pm
            );
    void evaluate_objects();
    void update_proj();
    void update();
private:
    int raw_feature_dim;
    InputList& inputs_;
    LabelList& labels_;
    ObjList& objects_;
    arma::mat& feature_base_;
    arma::mat& feature_centers_;
    std::vector<arma::uword> patch_frame_;
    std::vector<arma::uword> patch_class_;
    std::vector<double> patch_score_;
    arma::mat  patch_feature_;
    arma::mat  Sb;
    arma::mat  Sw;
    Config::Ptr config_;
};

#endif // UPDATECLUSTERCENTER_H