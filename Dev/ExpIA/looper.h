#ifndef LOOPER_H
#define LOOPER_H
#include <QObject>
#include "common.h"
#include "objectmodel.h"
#include "mainwindow.h"
class Looper:public QObject
{
    Q_OBJECT
public:
    Looper(
            MeshBundle<DefaultMesh>::PtrList& inputs,
            std::vector<ObjModel::Ptr>& objects,
            std::vector<arma::uvec>& labels,
            arma::mat& feature_base,
            arma::mat& feature_centers,
            MainWindow* main_window,
            QThread* main_thread,
            QObject* parent=0
           )
        :inputs_(inputs),
         objects_(objects),
         labels_(labels),
         feature_base_(feature_base),
         feature_centers_(feature_centers),
         main_window_(main_window),
         main_thread_(main_thread),
         QObject(parent){}
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
    void save_svx(QString);//save super voxel
    void save_lbl(QString);//save label
    void save_obj(QString);//save object model
    void save_clt(QString);//save cluster
    void showInMdi(QWidget* w, Qt::WindowFlags flag);
    void start_uo();
    void end();
public slots:
    void loop();
    void reset();
protected slots:
    void passMessage(QString msg,int t){emit message(msg,t);}
    void step_finished(){current_running_ = false;}
protected:
    void sv();//generate supervoxel
    void rg();//do region grow
    void uf();//unify label
    void uo();//update object model
    void uc();//upate cluster
    void gc();//graph cut
    void wait_for_current();
    void wait_for_current(QThread*);
protected:
    uint32_t count_;
    uint32_t max_count_;
    bool current_running_;
private:
    Config::Ptr config_;
    MainWindow* main_window_;
    QThread* main_thread_;
    MeshBundle<DefaultMesh>::PtrList& inputs_;
    std::vector<ObjModel::Ptr>& objects_;
    std::vector<arma::uvec>& labels_;
    arma::mat& feature_base_;
    arma::mat& feature_centers_;
};

#endif // LOOPER_H
