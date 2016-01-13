#ifndef UPDATECLUSTERCENTER_H
#define UPDATECLUSTERCENTER_H
#include <QThread>
#include "common.h"
class UpdateClusterCenter : public QThread
{
    Q_OBJECT
public:
    typedef std::shared_ptr<DefaultMesh> MeshPtr;
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> InputList;
    typedef std::vector<arma::uvec> LabelList;
    typedef std::vector<MeshPtr> ObjList;
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
              feature_centers_(center)
    {
        setObjectName("UpdateClusterCenter");
    }
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
protected:
    void run();
private:
    InputList& inputs_;
    LabelList& labels_;
    std::vector<MeshPtr>& objects_;
    arma::mat& feature_base_;
    arma::mat& feature_centers_;
    Config::Ptr config_;
};

#endif // UPDATECLUSTERCENTER_H
