#include "updateobjectmodel.h"
#include "ui_updateobjectmodel.h"
#include "registrationcore.h"
UpdateObjectModel::UpdateObjectModel(IMeshList &inputs, ILabelList &labels, OModelList &outputs, QWidget *parent) :
    QFrame(parent),
    inputs_(inputs),
    labels_(labels),
    outputs_(outputs),
    geo_thread_(NULL),
    color_thread_(NULL),
    ui(new Ui::UpdateObjectModel)
{
    ui->setupUi(this);
    geo_view_ = new MeshListViewerWidget();
    geo_view_->setMinimumSize(640,480);
    color_view_ = new LabSpace();
    ui->horizontalLayout->addWidget(geo_view_);
    ui->horizontalLayout->addWidget(color_view_);

    connect(&timer_,SIGNAL(timeout()),this,SLOT(prepare_for_next()));
    connect(&timer_,SIGNAL(timeout()),this,SLOT(start_align()));
    connect(&timer_,SIGNAL(timeout()),this,SLOT(start_fit()));

    connect(&gl_timer,SIGNAL(timeout()),geo_view_,SLOT(updateGL()));

    timer_.setSingleShot(true);
    gl_timer.setSingleShot(false);

    //find out the max label
    max_label_ = 0;
    ILabelList::iterator iter;
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        arma::uword max = arma::max(*iter);
        if(max>max_label_)max_label_=max;
    }
    current_label_ = 0;
    current_patches_.resize(inputs_.size());
    gl_timer.start(90);
}

bool UpdateObjectModel::configure(Config::Ptr config)
{
    config_ = config;
    if(!config_->has("Align_Max_Iter"))return false;
    if(!config_->has("Align_T_Eps"))return false;
    return true;
}

void UpdateObjectModel::startLater()
{
    timer_.start(10);
}

UpdateObjectModel::~UpdateObjectModel()
{
    geo_view_->close();
    color_view_->close();
    ui->horizontalLayout->removeWidget(geo_view_);
    ui->horizontalLayout->removeWidget(color_view_);
    geo_view_->deleteLater();
    color_view_->deleteLater();
    delete ui;
}

void UpdateObjectModel::prepare_for_next()
{
    ++current_label_;
    done_align_ = false;
    done_fit_ = false;
    extract_patches();
    if(valid_patches_.empty())
    {
        prepare_for_next();
    }
    geo_view_->list() = current_patches_;
}

void UpdateObjectModel::extract_patches()
{
    arma::uvec indices;
    ILabelList::iterator iter;
    IMeshList::iterator miter=inputs_.begin();
    PatchList::iterator piter=current_patches_.begin();
    QString tmp;
    size_t frame = 0;

    valid_patches_.clear();
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        indices = arma::find(*iter==current_label_);
        if(!*piter)*piter = std::make_shared<MeshBundle<DefaultMesh>>();
        tmp = tmp.sprintf("F%dL%d",frame,current_label_);
        (*piter)->name_ = tmp.toStdString();
        if(!indices.is_empty()){
            (*piter)->mesh_.clear();
            extractMesh<DefaultMesh,DefaultMesh>((*miter)->mesh_,indices,(*piter)->mesh_);
            valid_patches_.push_back(frame);
        }
        if(miter==inputs_.end())break;
        if(piter==current_patches_.end())break;
        ++miter;
        ++piter;
        ++frame;
    }
}

void UpdateObjectModel::update_objects()
{
    if(geo_thread_)
    {
        JRMPC_Thread* th = qobject_cast<JRMPC_Thread*>(geo_thread_);
        if(th)
        {
            Registration::JRMPC<DefaultMesh>::ResPtr r = th->result();
            while(outputs_.size()<current_label_)
            {
                outputs_.emplace_back(new ObjModel);
            }
            ObjModel::Ptr obj_ptr = outputs_[current_label_-1];
            if(!obj_ptr->GeoM_)obj_ptr->GeoM_.reset(new MeshBundle<DefaultMesh>);

            std::vector<arma::uword>::iterator iter;
            obj_ptr->GeoT_.resize(inputs_.size());
            //update transformation
            arma::uword index = 0;
            for(iter=valid_patches_.begin();iter!=valid_patches_.end();++iter)
            {
                ObjModel::T::Ptr& t_ptr = obj_ptr->GeoT_[*iter];
                t_ptr = std::make_shared<ObjModel::T>();
                arma::fmat R(t_ptr->R,3,3,false,true);
                arma::fvec t(t_ptr->t,3,false,true);
                R = *(r->Rs[index]);
                t = *(r->ts[index]);
                ++index;
            }
            //update geometry object
            std::vector<MeshBundle<DefaultMesh>::Ptr>& patch_list_ = geo_view_->list();
            PatchList::reverse_iterator piter = patch_list_.rbegin();
            ++piter;
            for(;piter!=patch_list_.rend();++piter)
            {
                MeshBundle<DefaultMesh>::Ptr& patch_ptr = *piter;
                if(!patch_ptr||0==patch_ptr.use_count())continue;
                obj_ptr->update(patch_ptr);
            }
            obj_ptr->computeLayout();;
        }
    }
}

void UpdateObjectModel::start_align()
{
    if(geo_thread_)
    {
        QString msg = "Please Wait Till the End of Last Registration\n";
        QMessageBox::critical(this, windowTitle(), msg);
        return;
    }
    JRMPC_Thread* th = new JRMPC_Thread();
    if(!th->init(geo_view_->list(),valid_patches_,config_))
    {
        QString msg = "Fail to Initialize the Registration:\n '";
        msg += QString::fromStdString(th->errorString());
        QMessageBox::critical( NULL, windowTitle(), msg);
        return;
    }else{
        geo_view_->reset_center();
        geo_view_->show_back();
        connect(th,SIGNAL(finished()),this,SLOT(finish_current()));
    }
    geo_thread_ = th;
    QString name;
    name = name.sprintf("Align L%d",current_label_);
    geo_thread_->setObjectName(name);
    emit message("Start "+name,2000);
    geo_thread_->start(QThread::HighestPriority);
}

void UpdateObjectModel::start_fit()
{
    ;
}

void UpdateObjectModel::finish_current()
{
    QThread* th = qobject_cast<QThread*>(sender());
    if(th)
    {
        QString msg = th->objectName() + " is Finished";
        if(th==geo_thread_)
        {
            while(geo_thread_->isRunning())
            {
                geo_thread_->terminate();
                QApplication::processEvents();
            }
            update_objects();
            show_layouts();
            done_align_ = true;
            geo_thread_->deleteLater();
            geo_thread_ = NULL;
        }
        if(th==color_thread_)
        {
            while(color_thread_->isRunning())
            {
                color_thread_->terminate();
                QApplication::processEvents();
            }
            done_fit_ = true;
            color_thread_->deleteLater();
            color_thread_ = NULL;
        }
        emit message(msg,0);
    }
    if(done_align_&&done_fit_)
    {
        timer_.start(1);
    }
}

void UpdateObjectModel::show_layouts(void)
{
    std::vector<ObjModel::T::Ptr>::iterator iter;
    ObjModel::Ptr obj_ptr = outputs_[current_label_-1];
    size_t index = 0;
    for(iter=obj_ptr->GeoT_.begin();iter!=obj_ptr->GeoT_.end();++iter)
    {
        ObjModel::T &T = **iter;
        arma::fmat::fixed<3,3> R(T.R);
        arma::fvec::fixed<3> t(T.t);
        arma::fmat invR = arma::inv(R);
        MeshBundle<DefaultMesh>::Ptr layout_m(new MeshBundle<DefaultMesh>);
        layout_m->mesh_ =
                    obj_ptr->GeoLayout_->mesh_;
        arma::fmat layout((float*)layout_m->mesh_.points(),3,layout_m->mesh_.n_vertices(),false,true);
        layout.each_col() -= t;
        layout = invR*layout;
        emit show_layout(index,layout_m);
        ++index;
    }
}
