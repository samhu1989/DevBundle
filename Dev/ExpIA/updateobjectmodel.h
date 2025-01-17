#ifndef UPDATEOBJECTMODEL_H
#define UPDATEOBJECTMODEL_H
#include <QFrame>
#include <armadillo>
#include <QTimer>
#include <QTime>
#include "common.h"
#include "objectmodel.h"
#include "MeshListViewerWidget.h"
#include "labspace.h"
#include "registrationcore.h"
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
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> PatchList;
    explicit UpdateObjectModel(
            IMeshList& inputs,
            ILabelList& labels,
            OModelList& outputs,
            QWidget *parent = 0);
    ~UpdateObjectModel();
    bool configure(Config::Ptr);
    void setMethod(uint32_t method){method_id_ = method;}
    inline void showMessageBox(bool show){show_mbox_=show;}
public slots:
    void startLater();
signals:
    void message(QString,int);
    void show_layout(int,MeshBundle<DefaultMesh>::Ptr);
    void closeInMdi(QWidget*);
protected slots:
    void prepare_for_next();
    void start_align(); //start thread for registration
    void start_fit();   //start thread for gmm fitting for color model
    void finish_current();
protected:
    JRMPC_Thread* create_align_thread();
    void extract_patches();
    void update_objects();
    void show_layouts();
private:
    Config::Ptr config_;
    IMeshList& inputs_;
    ILabelList& labels_;
    OModelList& outputs_;

    QTime alg_timer;
    QTimer timer_;
    QTimer gl_timer;
    uint32_t method_id_;
    QThread* geo_thread_;
    QThread* color_thread_;

    bool done_align_;
    bool done_fit_;
    bool show_mbox_;

    MeshListViewerWidget* geo_view_;
//    LabSpace* color_view_;
    Ui::UpdateObjectModel *ui;

    arma::uword current_label_;
    arma::uword max_label_;

    PatchList current_patches_;
    PatchList current_extanded_patches_;
    std::vector<arma::uword> valid_patches_;
};

#endif // UPDATEOBJECTMODEL_H
