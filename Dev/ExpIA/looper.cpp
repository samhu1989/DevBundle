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
bool Looper::configure(Config::Ptr config)
{
    config_= config;
    return true;
}

void Looper::reset()
{
    count_ = 0;
    max_count_ = 1;
}

void Looper::loop()
{
    QDir dir;
    reset();
    sv();
    saveStep(SV,0,0);
    for(count_=0;count_ < max_count_; ++count_ )
    {

        rg();
        saveStep(RG,count_,0);
        uf();
        saveStep(UF,count_,1);
        uo();
        saveStep(UO,count_,2);
        uc();
        saveStep(UC,count_,3);
        uf();
        saveStep(UF,count_,4);
        uo();
        saveStep(UO,count_,5);
        gc();
        saveStep(GC,count_,6);
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
    emit message(tr("Starting Update Object"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::gc()
{
    GraphCutThread* th = new GraphCutThread(inputs_,objects_,labels_);
    if(!th->configure(config_)){
        th->deleteLater();
        QString msg = "Missing Some Inputs or configure\n";
        QMessageBox::critical(NULL,tr("Looping:Failed to Start Graph Cut"), msg);
        return;
    }
    connect(th,SIGNAL(message(QString,int)),this,SLOT(passMessage(QString,int)));
    connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),main_window_,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
    emit message(tr("Starting Graph Cut"),0);
    th->start(QThread::HighPriority);
    wait_for_current(th);
    th->deleteLater();
}

void Looper::wait_for_current()
{
    while(current_running_)
    {
        QThread::msleep(30);
        QApplication::processEvents();
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
    }
    dir.mkdir(step_path);
    switch(step)
    {
    case SV:
        emit save_svx(dir.absoluteFilePath(step_path));break;
    case RG:
        emit save_lbl(dir.absoluteFilePath(step_path));break;
    case UF:
        emit save_lbl(dir.absoluteFilePath(step_path));break;
        emit save_clt(dir.absoluteFilePath(step_path));break;
    case UC:
        emit save_clt(dir.absoluteFilePath(step_path));break;
    case UO:
        emit save_obj(dir.absoluteFilePath(step_path));break;
    }
}

