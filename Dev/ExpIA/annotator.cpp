#include "annotator.h"
#include "ui_annotator.h"
#include "meshpairviewerwidget.h"
#include <QImage>
Annotator::Annotator(
        std::vector<QWidget*>&inputs,
        QTimer& gl_timer,
        std::vector<arma::uvec>&labels,
        QWidget *parent):
    inputs_(inputs),
    gl_timer_(gl_timer),
    labels_(labels),
    QFrame(parent),
    max_queue_size_(32),
    ui(new Ui::Annotator)
{
    ui->setupUi(this);
    connect(ui->label,SIGNAL(valueChanged(int)),this,SLOT(updateLabel(int)));
    connect(ui->apply2patch,SIGNAL(clicked(bool)),this,SLOT(applyToPatch()));
    connect(ui->apply2point,SIGNAL(clicked(bool)),this,SLOT(applyToPoint()));
    connect(ui->picklabel,SIGNAL(clicked(bool)),this,SLOT(pickLabel()));
}

bool Annotator::configure(Config::Ptr config)
{
    if(inputs_.empty())return false;
    if(labels_.empty())return false;
    if(inputs_.size()!=labels_.size())return false;
    init();
    return true;
}

void Annotator::init(void)
{
    uint64_t index = 0;
    for(std::vector<QWidget*>::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        arma::uvec& target = labels_[index];
        for(arma::uvec::iterator iter=target.begin();iter!=target.end();++iter)
        {
            updateLabel(int(*iter));
        }
        update_color(w->first_ptr(),target);
        ++index;
    }
}

void Annotator::update_color(MeshBundle<DefaultMesh>::Ptr m,const arma::uvec& label)
{
    assert(m->mesh_.n_vertices()==label.size());
    uint32_t* c = (uint32_t*)m->custom_color_.vertex_colors();
    #pragma omp parallel for
    for(int i=0;i<label.size();++i)
    {
        assert( l2c_.find(int(label(i))) != l2c_.end() );
        c[i] = l2c_[ int(label(i)) ];
    }
}

void Annotator::undo(void)
{
    undo_annotation();
}

void Annotator::pickLabel(void)
{
    uint64_t index = 0;
    bool picked = false;
    int label = -1;
    for(std::vector<QWidget*>::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        if(picked)w->first_selected().clear();
        else {
            if(!w->first_selected().empty())
            {
                label = labels_[index](w->first_selected().back());
                picked = true;
                w->first_selected().clear();
            }
        }
        QApplication::processEvents();
        ++index;
    }
    if( label >= 0 )
    {
        updateLabel(label);
    }
}

void Annotator::applyToPatch(void)
{
//    std::cerr<<"applying to patch"<<std::endl;
    uint64_t index = 0;
    int label = ui->label->value();
    for(std::vector<QWidget*>::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        if(!w->first_selected().empty())
        {
            arma::uvec& target = labels_[index];
            std::vector<arma::uword>::iterator e = w->first_selected().end();
            std::vector<arma::uword> selected_patch;
            selected_patch.reserve(w->first_selected().size());
            for(std::vector<arma::uword>::iterator iter=w->first_selected().begin();iter!=e;++iter)
            {
                selected_patch.push_back( target(*iter) );
            }
            std::sort(selected_patch.begin(),selected_patch.end());
            selected_patch.erase(std::unique(selected_patch.begin(),selected_patch.end()),selected_patch.end());
            arma::uvec indices = arma::find(target==selected_patch.front());
            for(std::vector<arma::uword>::iterator iter = selected_patch.begin()+1 ; iter!=selected_patch.end() ; ++iter )
            {
                indices = arma::join_rows(indices,arma::find(target==*iter));
            }
            w->first_selected().clear();
            annotation_deque_.emplace_back(new Annotation(indices,label,w->first_ptr(),target));
            annotation_deque_.back()->apply();
            if( annotation_deque_.size() > max_queue_size_ )
            {
                annotation_deque_.pop_front();
            }
            update_color(w->first_ptr(),target);
            QApplication::processEvents();
        }
        ++index;
    }
}

void Annotator::applyToPoint(void)
{
    uint64_t index = 0;
    int label = ui->label->value();
    for(std::vector<QWidget*>::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        if(!w->first_selected().empty())
        {
            arma::uvec& target = labels_[index];
            arma::uvec indices(w->first_selected());
            w->first_selected().clear();
            annotation_deque_.emplace_back(new Annotation(indices,label,w->first_ptr(),target));
            annotation_deque_.back()->apply();
            if( annotation_deque_.size() > max_queue_size_ )
            {
                annotation_deque_.pop_front();
            }
            update_color(w->first_ptr(),target);
            QApplication::processEvents();
        }
        ++index;
    }
}

void Annotator::updateLabel(int label)
{
    if( ui->label->value() != label )
    {
        ui->label->blockSignals(true);
        ui->label->setValue(label);
        ui->label->blockSignals(false);
    }
    QRgb c;
    if(l2c_.find(ui->label->value())!=l2c_.end())
    {
        c = l2c_[ui->label->value()];
    }else{
        c = ColorArray::rand_color();
        while(c2l_.find(c)!=c2l_.end())
        {
            c = ColorArray::rand_color();
        }
        c2l_[c] = ui->label->value();
        l2c_[ui->label->value()] = c;
    }
    QImage img(32,32,QImage::Format_RGB888);
    QColor qc(c);
    img.fill(QColor(qc.blue(),qc.green(),qc.red()));
    ui->labelcolor->setPixmap(QPixmap::fromImage(img));
}

void Annotator::undo_annotation()
{
    if(!annotation_deque_.empty())
    {
        annotation_deque_.back()->undo();
        update_color(annotation_deque_.back()->mesh_ptr(),annotation_deque_.back()->label());
        annotation_deque_.pop_back();
    }else {
        emit message(tr("No annotation in queue"),3000);
    }
}

Annotator::~Annotator()
{
    delete ui;
}
