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
    reset();
    sv();
    for(count_=0;count_ < max_count_; ++count_ )
    {
        rg();
        uf();
        uo();
        uc();
        uf();
        uo();
        gc();
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
    wait_for_current(th);
    th->deleteLater();
}

void Looper::uo()
{
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
    wait_for_current(th);
    th->deleteLater();
}

void Looper::wait_for_current()
{
    while(current_running_)
    {
        QThread::msleep(300);
        QApplication::processEvents();
    }
}

void Looper::wait_for_current(QThread* th)
{
    while(!th->isFinished())
    {
        QThread::msleep(300);
        QApplication::processEvents();
    }
}

