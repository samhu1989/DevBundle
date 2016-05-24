#ifndef JRCSVIEW_H
#define JRCSVIEW_H

#include <QFrame>
#include "common.h"
#include "objectmodel.h"
#include "MeshListViewerWidget.h"
#include "jrcsthread.h"
#include <QTimer>
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
private:
    Ui::JRCSView *ui;
    MeshListViewerWidget* geo_view_;
    QThread* jrcs_thread_;
    JRCSThread* jrcs_worker_;
    JRCS::JRCSBase::TsLst rt_;
    MeshList& inputs_;
    LabelList& labels_;
    ModelList& objects_;
    QTimer t_;
};

#endif // JRCSVIEW_H
