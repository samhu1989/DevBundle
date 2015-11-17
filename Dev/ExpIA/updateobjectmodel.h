#ifndef UPDATEOBJECTMODEL_H
#define UPDATEOBJECTMODEL_H
#include <QFrame>
#include <armadillo>
#include "common.h"
#include "objectmodel.h"
namespace Ui {
class UpdateObjectModel;
}

class UpdateObjectModel : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> IMeshList;
    typedef std::vector<arma::uvec> ILabelList;
    typedef std::vector<ObjModel::Ptr> OModelList;
    explicit UpdateObjectModel(
            IMeshList& inputs,
            ILabelList& labels,
            OModelList& outputs,
            QWidget *parent = 0);
    ~UpdateObjectModel();
    bool configure(Config::Ptr);
public slots:
    void startLater();
signals:
    void message(QString,int);
protected slots:

private:
    Config::Ptr config_;
    IMeshList& inputs_;
    ILabelList& labels_;
    OModelList& outputs_;
    Ui::UpdateObjectModel *ui;
};

#endif // UPDATEOBJECTMODEL_H
