#ifndef JRCSVIEW_H
#define JRCSVIEW_H

#include <QFrame>
#include "common.h"
#include "objectmodel.h"
#include "MeshListViewerWidget.h"
#include "jrcsthread.h"
#include <QTimer>
#include <QTime>
namespace Ui {
class JRCSView;
}

class JRCSView : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    typedef std::vector<arma::uvec> LabelList;
    typedef std::vector<ObjModel::Ptr> ModelList;
    typedef JRCSThread::MatPtr MatPtr;
    typedef JRCSThread::CMatPtr CMatPtr;
    typedef JRCSThread::LMatPtr LMatPtr;
    typedef JRCSThread::MatPtrLst MatPtrLst;
    typedef JRCSThread::CMatPtrLst CMatPtrLst;
    typedef JRCSThread::LCMatPtrLst LCMatPtrLst;
    typedef JRCSThread::LMatPtrLst LMatPtrLst;
    explicit JRCSView(
            MeshList& inputs,
            LabelList& labels,
            ModelList& objects,
            QWidget *parent = 0
            );
    ~JRCSView();
    bool configure(Config::Ptr);
    bool init(Config::Ptr);
    inline void set_method(std::shared_ptr<JRCS::JRCSBase> method){jrcs_worker_->set_method(method);}
    inline void set_init_method(std::shared_ptr<JRCS::JRCSInitBase> method){jrcs_worker_->set_init_method(method);}
    inline void set_show_mode(const std::string& _s){geo_view_->set_draw_mode(_s);}
    void start();
signals:
    void message(QString,int);
    void closeInMdi(QWidget*);
protected slots:
    void passMessage(QString,int);
    void finished();
    void align(int);
protected:
    void input( JRCSThread* jrcs_worker_ );
    bool allocate_x( JRCSThread* jrcs_worker_ );
    void move_worker_to_thread( JRCSThread* jrcs_worker );
    void save_rt();
    void save_centroids();
    void save_order();//corrrespondence by order
private:
    Ui::JRCSView *ui;
    MeshListViewerWidget* geo_view_;
    QThread* jrcs_thread_;
    JRCSThread* jrcs_worker_;
    JRCS::JRCSBase::TsLst rt_;
    std::vector<arma::uvec> order_;
    MeshList& inputs_;
    LabelList& labels_;
    ModelList& objects_;
    QTimer t_;
    QTime time;
    bool close_on_finish_;
};

#endif // JRCSVIEW_H
