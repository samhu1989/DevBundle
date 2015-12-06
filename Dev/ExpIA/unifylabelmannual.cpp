#include "unifylabelmannual.h"
#include "ui_unifylabelmannual.h"
#include <memory>
#include <QScrollBar>
UnifyLabelMannual::UnifyLabelMannual(
        InputType&input,
        OutputType&output,
        QWidget*parent
        ) :
    QFrame(parent),
    current_frame_(0),
    current_patch_(0),
    inputs_(input),
    labels_(output),
    ui(new Ui::UnifyLabelMannual)
{
    ui->setupUi(this);
    reloadFrame();
    connect(&timer,SIGNAL(timeout()),this,SLOT(init()));
    timer.setSingleShot(true);
}

bool UnifyLabelMannual::configure(Config::Ptr config)
{
    config_ = config;
    OutputType::iterator iter;
    for( iter = labels_.begin() ; iter != labels_.end() ; ++iter )
    {
        if(0==arma::max(*iter)||iter->is_empty()){
            std::clog<<"You probably should do regiongrow first\n";
            return false;
        }
    }
    return true;
}

void UnifyLabelMannual::init(void)
{
    showPatches();
}

void UnifyLabelMannual::keyPressEvent(QKeyEvent* e)
{
    std::cerr<<e->text().toStdString()<<std::endl;
    switch(e->key())
    {
    case Qt::Key_W:
        frameLast();
        break;
    case Qt::Key_S:
        frameNext();
        break;
    case Qt::Key_A:
        movePatchLast();
        break;
    case Qt::Key_D:
        movePatchNext();
        break;
    case Qt::Key_Space:
        patchNext();
        break;
    case Qt::Key_Delete:
        patchDelete();
        break;
    default:
        showHelp();
    }
    QFrame::keyPressEvent(e);
}

void UnifyLabelMannual::frameLast()
{
    updateObjects();
    updateLabel();
    if( current_frame_ > 0 )current_frame_--;
    else current_frame_ = inputs_.size() - 1;
    reloadFrame();
    showPatches();
}

void UnifyLabelMannual::frameNext()
{
    updateObjects();
    updateLabel();
    if( current_frame_ < inputs_.size() - 1 )current_frame_++;
    else current_frame_ = 0;
    reloadFrame();
    showPatches();
}

void UnifyLabelMannual::movePatchNext()
{
    size_t next_patch_;
    if( current_patch_ == patches_.size()-1 )
    {
        //insertion is allowed to the end if the last object is filled
        PatchPairView* v = pair_views_[current_patch_];
        if( v->oview()->mesh_ptr() && ( 0 != v->oview()->mesh_ptr()->n_vertices() ) )
        {
            patches_.push_back(std::make_shared<DefaultMesh>());
            patch_label_.insert_rows(patch_label_.size(),1);
            patch_label_(patch_label_.size()-1) = patches_.size();
            PatchPairView* w = new PatchPairView();
            w->setAttribute(Qt::WA_DeleteOnClose,true);
            pair_views_.push_back(w);
            ui->layout->addWidget(w);
            pair_views_.back()->show();
            next_patch_ = current_patch_ + 1;
        }else return;
    }
    else next_patch_ = current_patch_ + 1;
    //exchange in patch
    MeshPtr ptr = patches_[current_patch_];
    patches_[current_patch_] = patches_[next_patch_];
    patches_[next_patch_] = ptr;
    //exchange in label
    arma::uword l = patch_label_[current_patch_];
    patch_label_[current_patch_] = patch_label_[next_patch_];
    patch_label_[next_patch_] = l;
    //reset in view
    pair_views_[current_patch_]->pview()->reset_ptr(patches_[current_patch_]);
    pair_views_[next_patch_]->pview()->reset_ptr(patches_[next_patch_]);
    //update selection
    pair_views_[current_patch_]->setFrameShape(QFrame::NoFrame);
    pair_views_[next_patch_]->setFrameShape(QFrame::Box);
    //
    current_patch_ = next_patch_;

}

void UnifyLabelMannual::movePatchLast()
{
    size_t last_patch_;
    if(current_patch_ == 0)return;
    else last_patch_ = current_patch_ - 1;
    //exchange in patch
    MeshPtr ptr = patches_[current_patch_];
    patches_[current_patch_] = patches_[last_patch_];
    patches_[last_patch_] = ptr;
    //exchange in label
    arma::uword l = patch_label_[current_patch_];
    patch_label_[current_patch_] = patch_label_[last_patch_];
    patch_label_[last_patch_] = l;
    //reset in view
    pair_views_[current_patch_]->pview()->reset_ptr(patches_[current_patch_]);
    pair_views_[last_patch_]->pview()->reset_ptr(patches_[last_patch_]);
    //update selection
    pair_views_[current_patch_]->setFrameShape(QFrame::NoFrame);
    pair_views_[last_patch_]->setFrameShape(QFrame::Box);
    //
    current_patch_ = last_patch_;
}

void UnifyLabelMannual::patchNext()
{
    PatchPairView* v = pair_views_[current_patch_];
    v->setFrameShape(QFrame::NoFrame);
    if( current_patch_ < patches_.size() - 1 )++current_patch_;
    else current_patch_ = 0;
    v = pair_views_[current_patch_];
    v->setFrameShape(QFrame::Box);
    QScrollBar* bar = ui->scrollArea->horizontalScrollBar();
    int value;
    int min = bar->minimum();
    int max = bar->maximum();
    value = current_patch_ * ( max - min ) / ( patches_.size() - 1 );
    bar->setValue(value+min);
}

void UnifyLabelMannual::patchDelete()
{
    PatchPairView* v = pair_views_[current_patch_];
    if(v->oview()->mesh_ptr()&&(0!=v->oview()->mesh_ptr()->n_vertices()))//empty view
    {
        MeshPtr ptr(new DefaultMesh);
        v->pview()->reset_ptr(ptr);
        patches_[current_patch_]->clear();
        patch_label_[current_patch_] = 0;
    }else{
        if(current_patch_==patches_.size()-1)
        {
            patches_.pop_back();
            patch_label_.shed_row(patch_label_.size()-1);
            --current_patch_;
        }else
        {
            patches_.erase(patches_.begin()+current_patch_);
            patch_label_.shed_row(current_patch_);
            for(int i=current_patch_;i<pair_views_.size()-1;++i)
            {
                pair_views_[i]->pview()->reset_ptr(pair_views_[i+1]->pview()->mesh_ptr());
            }
        }
        pair_views_.back()->hide();
        pair_views_.back()->close();
        ui->layout->removeWidget(pair_views_.back());
        pair_views_.back()->deleteLater();
        update();
        pair_views_.pop_back();
        pair_views_[current_patch_]->setFrameShape(QFrame::Box);
    }
}

void UnifyLabelMannual::showHelp()
{
    ;
}

void UnifyLabelMannual::reloadFrame()
{
    QString msg;
    emit message(msg.sprintf("Current Frame: %d",current_frame_),0);
    MeshBundle<DefaultMesh>& bundle = *inputs_[current_frame_];
    arma::uvec& label = labels_[current_frame_];
    arma::uword maxL = arma::max(label);
    arma::uword l;
    arma::uvec indices;
    patches_.clear();
    l = 1;
    while( (patches_.size() < pair_views_.size()) || (l <= maxL) )
    {
        patches_.push_back(std::make_shared<DefaultMesh>());
        if( l <= maxL )
        {
            indices = arma::find( label == l );
            if(!indices.is_empty())extractMesh<DefaultMesh,DefaultMesh>(bundle.mesh_,indices,*patches_.back());
        }//else let the patch be empty mesh
        l++;
    }
    if( current_patch_ >= patches_.size() )
    {
        if(!pair_views_.empty())pair_views_[current_patch_]->setFrameShape(QFrame::NoFrame);
        current_patch_ = 0;
    }
    patch_label_ = arma::linspace<arma::uvec>(1,patches_.size(),patches_.size());
}

void UnifyLabelMannual::showPatches()
{
    while( pair_views_.size() < patches_.size() )
    {
        PatchPairView* w = new PatchPairView();
        w->setAttribute(Qt::WA_DeleteOnClose,true);
        pair_views_.push_back(w);
        ui->layout->addWidget(w);
        pair_views_.back()->show();
    }
    std::vector<PatchPairView*>::iterator witer;
    std::vector<MeshPtr>::iterator miter = patches_.begin();
    for(witer=pair_views_.begin();witer!=pair_views_.end();++witer)
    {
        (*witer)->pview()->reset_ptr(*miter);
        ++miter;
        if(miter==patches_.end())break;
    }
    pair_views_[current_patch_]->setFrameShape(QFrame::Box);
}

void UnifyLabelMannual::updateLabel()
{
    arma::uvec::iterator iter;
    arma::uvec& current_label_ = labels_[current_frame_];
    arma::uvec new_label_(current_label_.size(),arma::fill::zeros);
    arma::uword label = 1;
    for( iter = patch_label_.begin() ;iter != patch_label_.end() ; ++iter )
    {
        if(*iter!=0)//a patch label start from 1
        {
            new_label_.elem(arma::find(current_label_==*iter)).fill(label);
        }
        ++label;
    }
    current_label_ = new_label_;
    inputs_[current_frame_]->custom_color_.fromlabel(current_label_);
}

void UnifyLabelMannual::updateObjects()
{
    while( objects_.size() < pair_views_.size() )
    {
        PatchPairView* v = pair_views_[objects_.size()];
        objects_.push_back(patches_[objects_.size()]);
        v->oview()->reset_ptr(objects_.back());
    }
}

UnifyLabelMannual::~UnifyLabelMannual()
{
    delete ui;
}
