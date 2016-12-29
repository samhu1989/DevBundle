#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "updateclustercenter.h"
#include "looper.h"
#include "jrcsview.h"
#include "jrcsinitthread.h"
#include "regiongrowrgbthread.h"
#include "sort_agd.h"
#include "ncut.h"
#include "supervoxelthread.h"
#include "graphcutthread.h"
#include "regiongrowthread.h"
#include "unifylabelthread.h"
#include "unifylabelmannual.h"
#include "updateobjectmodel.h"
#include "globalalign.h"
#include "extractbackground.h"
#include "inpatchgraphcut.h"
#include "robustcut.h"
#include "jrcs.h"
#include "labelcompactor.h"
#include "annotator.h"
#include "gdcthread.h"
#include "jrcsprimitive.h"
void MainWindow::start_editing()
{
    if(edit_thread_)
    {
        QString msg = "Please Wait Till the End of Last Algorithm\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    if( inputs_.empty() || labels_.empty() )
    {
        QString msg = "No Inputs for Editing\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    QAction* edit = qobject_cast<QAction*>(sender());
    if(edit==ui->actionSupervoxel)
    {
        SupervoxelThread* th = new SupervoxelThread(inputs_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionRegionGrow)
    {
        RegionGrowThread* th = new RegionGrowThread(inputs_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionAutomatically)
    {
        UnifyLabelThread* th = new UnifyLabelThread(inputs_,labels_,feature_base_,feature_centers_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionMannually)
    {
        UnifyLabelMannual* w = new UnifyLabelMannual(inputs_,labels_);
        if(!w->configure(config_)){
            QString msg = "You probably should do regiongrow first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        showInMdi((QWidget*)w);
        w->initLater();
    }
    if(edit==ui->actionV0)
    {
        UpdateObjectModel* w = new UpdateObjectModel(
                    inputs_,
                    labels_,
                    objects_
                    );
        if(!w->configure(config_)){
            QString msg = "You probably should do unify label first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        w->setMethod(0);
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(w,SIGNAL(show_layout(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,Qt::Widget|Qt::WindowMinMaxButtonsHint);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->startLater();
    }
    if(edit==ui->actionV1)
    {
        UpdateObjectModel* w = new UpdateObjectModel(
                    inputs_,
                    labels_,
                    objects_
                    );
        if(!w->configure(config_)){
            QString msg = "You probably should do unify label first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        w->setMethod(1);
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(w,SIGNAL(show_layout(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,Qt::Widget|Qt::WindowMinMaxButtonsHint);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->startLater();
    }
    if(edit==ui->actionUpdate_Cluster_Center)
    {
        UpdateClusterCenter* th = new UpdateClusterCenter(inputs_,labels_,objects_,feature_base_,feature_centers_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        edit_thread_ = th;
    }
    if(edit==ui->actionGlobal_Align)
    {
        GlobalAlign* w = new GlobalAlign(inputs_);
        if(!w->configure(config_)){
            QString msg = "You probably should open inputs first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        s->show();
    }
    if(edit==ui->actionExtract_Background)
    {
        ExtractBackground* w = new ExtractBackground(mesh_views_,gl_timer,labels_);
        if(!w->configure(config_)){
            QString msg = "You probably should open inputs first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        s->show();
    }
    if(edit==ui->actionAnnotator)
    {
        Annotator* w = new Annotator(mesh_views_,gl_timer,labels_);
        if(!w->configure(config_)){
            QString msg = "You probably should open inputs first\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(ui->actionUndo,SIGNAL(triggered(bool)),w,SLOT(undo()));
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        s->show();
    }
    if(edit==ui->actionGlobal_Graph_Cut)
    {
        GraphCutThread* th = new GraphCutThread(inputs_,objects_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
        edit_thread_ = th;
    }
    if(edit==ui->actionIn_Patch_Graph_Cut)
    {
        InPatchGraphCut* th = new InPatchGraphCut(inputs_,objects_,labels_);
        if(!th->configure(config_)){
            th->deleteLater();
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        connect(th,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(th,SIGNAL(sendMatch(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
        edit_thread_ = th;
    }
    if(edit==ui->actionRegionGrowRGB)
    {
        RegionGrowRGBThread* worker = new RegionGrowRGBThread(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(finished()),th,SLOT(quit()));
        connect(worker,SIGNAL(finished()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("RegionGrowRGBThread"));
        edit_thread_ = th;
    }
    if(edit==ui->actionCut_Graph)
    {
        NCut* worker = new NCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("NCut"));
        edit_thread_ = th;
    }
    if(edit==ui->actionGet_Compact_Label)
    {
        LabelCompactor* worker = new LabelCompactor(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("Get Compact Label"));
        edit_thread_ = th;
    }
    if(edit==ui->actionDebug_Convexity)
    {
        NCut* worker = new NCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(debug_convexity()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("NCut_Debug_Convexity"));
        edit_thread_ = th;
    }
    if(edit==ui->actionDebug_Color)
    {
        NCut* worker = new NCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(debug_color()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("NCut_Debug_Color"));
        edit_thread_ = th;
    }
    if(edit==ui->actionDebug_Dist)
    {
        NCut* worker = new NCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(debug_dist()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("NCut_Debug_Dist"));
        edit_thread_ = th;
    }
    if(edit==ui->actionDebug_W)
    {
        NCut* worker = new NCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(debug_W()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("NCut_Debug_W"));
        edit_thread_ = th;
    }
    if(edit==ui->actionBase_Segments)
    {
        RobustCut* worker = new RobustCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(base_segments()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(this,SIGNAL(keyPressSignal(QKeyEvent*)),worker,SLOT(keyPressEvent(QKeyEvent*)));
        th->setObjectName(tr("Base_Segments"));
        edit_thread_ = th;
    }
    if(edit==ui->actionConsensus_Segment)
    {
        RobustCut* worker = new RobustCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(consensus_segment()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(this,SIGNAL(keyPressSignal(QKeyEvent*)),worker,SLOT(keyPressEvent(QKeyEvent*)));
        th->setObjectName(tr("Consensus_Segment"));
        edit_thread_ = th;
    }
    if(edit==ui->actionShow_Base_Segments)
    {
        RobustCut* worker = new RobustCut(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(show_base_segment()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(this,SIGNAL(keyPressSignal(QKeyEvent*)),worker,SLOT(keyPressEvent(QKeyEvent*)));
        th->setObjectName(tr("Show_Base_Segments"));
        edit_thread_ = th;
    }
    if(edit==ui->actionSort_AGD)
    {
        Sort_AGD* worker = new Sort_AGD(inputs_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(finished()),th,SLOT(quit()));
        connect(worker,SIGNAL(finished()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("Sort_AGD"));
        edit_thread_ = th;
    }
    if(edit==ui->actionJRCS_Init_SIHKS)
    {
        JRCSWork* worker = new JRCSWork(inputs_,labels_,objects_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(Init_SI_HKS()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("JRCS_Init_SI_HKS"));
        edit_thread_ = th;
    }
    if(edit==ui->actionDebug_HKS_Clustering)
    {
        JRCSWork* worker = new JRCSWork(inputs_,labels_,objects_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(debug_SI_HKS()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("JRCS_Debug_SI_HKS"));
        edit_thread_ = th;
    }
    if(edit==ui->actionJRCS_Init_Bernoulli)
    {
        JRCSWork* worker = new JRCSWork(inputs_,labels_,objects_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(Init_Bernolli()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("JRCS_Init_Bernoulli"));
        edit_thread_ = th;
    }
    if(edit==ui->actionJRCS_Init)
    {
        JRCSInitThread* worker = new JRCSInitThread(inputs_,labels_);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(finished()),th,SLOT(quit()));
        connect(worker,SIGNAL(finished()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(worker,SIGNAL(showbox(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)),Qt::DirectConnection);
        th->setObjectName(tr("JRCS_Init"));
        edit_thread_ = th;
    }
    if(edit==ui->actionJRCS_Opt_Basic)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        if(!JRCSWork::init_optimize(w))
        {
            QString msg = "Please run one of the external initializations\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionJRCS_Opt_AONI)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        JRCSWork::set_opt_aoni(w);
        if(!JRCSWork::init_optimize(w))
        {
            QString msg = "Please run one of the external initializations\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionJRCS_Opt_AOPT)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        JRCSWork::set_opt_aopt(w);
        if(!JRCSWork::init_optimize(w))
        {
            QString msg = "Please run one of the external initializations\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionJRCS_Opt_Spectrum)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        JRCSWork::set_opt_spectrum(w);
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionJRCS_Opt_Bilateral)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        JRCSWork::set_opt_bilateral(w);
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionJRCS_Opt_Primitive)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        std::shared_ptr<JRCS::JRCSBase> method(new JRCS::JRCSPrimitive());
        w->set_method(method);
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->set_show_mode("Colored Vertices");
//        w->start();
    }
    if(edit==ui->actionJRCS_Old)
    {
        JRCSView* w = new JRCSView(
                    inputs_,
                    labels_,
                    objects_
                    );
        if(!w->configure(config_)){
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            w->deleteLater();
            return;
        }
        connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        QMdiSubWindow* s = ui->mdiArea->addSubWindow(w);
        connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
        s->show();
        w->start();
    }
    if(edit==ui->actionIterate)
    {

        QThread* th = new QThread();
        Looper* loop = new Looper(inputs_,objects_,labels_,feature_base_,feature_centers_,this,QThread::currentThread());
        if(!loop->configure(config_)){
            loop->deleteLater();
            QString msg = "Missing Some Inputs or configure Before Iterate\n";
            QMessageBox::critical(this, windowTitle(), msg);
            return;
        }
        QString dirName = QFileDialog::getExistingDirectory(
                this,
                tr("Save Steps"),
                tr("../Dev_Data/")
                );
        if(dirName.isEmpty())return;
        loop->setSavePath(dirName);
        loop->moveToThread(th);
        connect(loop,SIGNAL(save_svx(QString)),this,SLOT(save_supervoxels(QString)));
        connect(loop,SIGNAL(save_obj(QString)),this,SLOT(save_objects(QString)));
        connect(loop,SIGNAL(save_lbl(QString)),this,SLOT(save_labels(QString)));
        connect(loop,SIGNAL(save_clt(QString)),this,SLOT(save_cluster(QString)));
        connect(loop,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        connect(loop,SIGNAL(end()),th,SLOT(quit()));
        connect(loop,SIGNAL(end()),loop,SLOT(deleteLater()));
        connect(th,SIGNAL(started()),loop,SLOT(loop()));
        connect(this,SIGNAL(object_updated()),loop,SLOT(step_finished()));
        connect(loop,SIGNAL(start_uo()),this,SLOT(update_object()));
        th->setObjectName(QString("Global_Iterate"));

        edit_thread_ = th;
    }

    if(edit==ui->actionGDCoord)
    {
        GDCThread* worker = new GDCThread(mesh_views_,gl_timer);
        if(!worker->configure(config_))
        {
            QString msg = "Missing Some Inputs or configure\n";
            QMessageBox::critical(this, windowTitle(), msg);
            worker->deleteLater();
            return;
        }
        QThread* th = new QThread();
        worker->moveToThread(th);
        connect(th,SIGNAL(started()),worker,SLOT(process()));
        connect(worker,SIGNAL(end()),th,SLOT(quit()));
        connect(worker,SIGNAL(end()),worker,SLOT(deleteLater()));
        connect(worker,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
        th->setObjectName(tr("GDCoord"));
        edit_thread_ = th;
    }

    if(edit_thread_){
        connect(edit_thread_,SIGNAL(finished()),this,SLOT(finish_editing()));
        edit_thread_->start(QThread::HighestPriority);
    }
}

void MainWindow::update_object()
{
    UpdateObjectModel* w = new UpdateObjectModel(
                inputs_,
                labels_,
                objects_
                );
    if(!w->configure(config_)){
        QString msg = "You probably should do unify label first\n";
        QMessageBox::critical(this, windowTitle(), msg);
        w->deleteLater();
        return;
    }
    w->setMethod(1);
    w->showMessageBox(false);
    connect(w,SIGNAL(message(QString,int)),ui->statusBar,SLOT(showMessage(QString,int)));
    connect(w,SIGNAL(show_layout(int,MeshBundle<DefaultMesh>::Ptr)),this,SLOT(showBox(int,MeshBundle<DefaultMesh>::Ptr)));
    w->setAttribute(Qt::WA_DeleteOnClose,true);
    QMdiSubWindow* s = ui->mdiArea->addSubWindow(w,Qt::Widget|Qt::WindowMinMaxButtonsHint);
    connect(w,SIGNAL(closeInMdi(QWidget*)),this,SLOT(closeInMdi(QWidget*)));
    s->show();
    connect(w,SIGNAL(destroyed(QObject*)),this,SLOT(notify_object_updated()));
    w->startLater();
}

void MainWindow::finish_editing()
{
    QString msg = edit_thread_->objectName() + " is Finished";
    if(edit_thread_)
    {
        while(edit_thread_->isRunning())
        {
            edit_thread_->quit();
            QApplication::processEvents();
        }
        edit_thread_->deleteLater();
        edit_thread_ = NULL;
    }
    QMessageBox::information( this, windowTitle(), msg);
}
