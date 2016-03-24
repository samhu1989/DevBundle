#include "looper.h"
#include <QDebug>
#include "updateobjectmodel.h"
#include <QMessageBox>
#include <QApplication>
#include "supervoxelthread.h"
#include "regiongrowthread.h"
#include "unifylabelthread.h"
#include "graphcutthread.h"
#include "updateclustercenter.h"
#include "inpatchgraphcut.h"
bool Looper::configure(Config::Ptr config)
{
    config_= config;
    if(config_->has("Loop_iter_num"))
    {
        max_count_ = config->getInt("Loop_iter_num");
    }else return false;
    return true;
}

void Looper::reset()
{
    count_ = 0;
    max_count_ = 1;
}

void Looper::loop()
{
    reset();
    sv();
    saveStep(SV,0,0);
    for(count_=0;count_ < max_count_; ++count_ )
    {
        rg();
        saveStep(RG,count_,1);
        uf();
        saveStep(UF,count_,2);
        uo();
        saveStep(UO,count_,3);
        uc();
        saveStep(UC,count_,4);
        uf();
        saveStep(UF,count_,5);
        uo();
        saveStep(UO,count_,6);
        ggc();
        saveStep(GGC,count_,7);
    }
    emit end();
}

void Looper::sv()
{
    SupervoxelThread* th = new SupervoxelThread(inputs_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Supervoxel"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    emit message(tr("Starting Supervoxel"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::rg()
{
    RegionGrowThread* th = new RegionGrowThread(inputs_,labels_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Region Grow"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    emit message(tr("Starting Region Grow"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::uf()
{
    UnifyLabelThread* th = new UnifyLabelThread(inputs_,labels_,feature_base_,feature_centers_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Unify Label"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    emit message(tr("Starting Unify Label"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::uo()
{
    emit message(tr("Starting Update Object"),0);
    emit start_uo();
    wait_for_current();
}

void Looper::uc()
{
    UpdateClusterCenter* th = new UpdateClusterCenter(inputs_,labels_,objects_,feature_base_,feature_centers_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Inputs or configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Update Cluster"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    emit message(tr("Starting Update Cluster Center"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::ggc()
{
    GraphCutThread* th = new GraphCutThread(inputs_,objects_,labels_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Inputs or configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Global Graph Cut"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),main_window_,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
    emit message(tr("Starting Global Graph Cut"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::lgc()
{
    InPatchGraphCut* th = new InPatchGraphCut(inputs_,objects_,labels_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Inputs or configure\n";
        QMessageBox::critical(NULL, tr("Looping:Failed to Start In-Patch Graph Cut"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
    emit message(tr("Starting In-Patch Graph Cut"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::wait_for_current()
{
    current_running_ = true;
    while(current_running_)
    {
        QThread::msleep(30);
        QApplication::processEvents();//it is suppose to recieve a signal of finish
    }
}

void Looper::wait_for_current(QThread* th)
{
    while(!th->isFinished())
    {
        QThread::msleep(30);
        QApplication::processEvents();
    }
}

void Looper::saveStep(Step step,int iter,int stepn)
{
    QDir dir;
    dir.setPath(save_path_);
    QString step_path;
    switch(step)
    {
    case SV:
         step_path = step_path.sprintf("./%02u%02usv",iter,stepn);break;
    case RG:
         step_path = step_path.sprintf("./%02u%02urg",iter,stepn);break;
    case UF:
         step_path = step_path.sprintf("./%02u%02uuf",iter,stepn);break;
    case UC:
         step_path = step_path.sprintf("./%02u%02uuc",iter,stepn);break;
    case UO:
         step_path = step_path.sprintf("./%02u%02uuo",iter,stepn);break;
    case GGC:
         step_path = step_path.sprintf("./%02u%02uggc",iter,stepn);break;
    }
    dir.mkdir(step_path);
    switch(step)
    {
    case SV:
        emit save_svx(dir.absoluteFilePath(step_path));break;
    case RG:
        emit save_lbl(dir.absoluteFilePath(step_path));break;
    case UF:
        emit save_lbl(dir.absoluteFilePath(step_path));
        emit save_clt(dir.absoluteFilePath(step_path));
        break;
    case UC:
        emit save_clt(dir.absoluteFilePath(step_path));break;
    case UO:
        emit save_obj(dir.absoluteFilePath(step_path));break;
    case GGC:
        emit save_lbl(dir.absoluteFilePath(step_path));break;
    }
}

