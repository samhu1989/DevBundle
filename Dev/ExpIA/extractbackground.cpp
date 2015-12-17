#include "extractbackground.h"
#include "ui_extractbackground.h"
#include "MeshPairViewerWidget.h"
#include "segmentationcore.h"
#include "extractplanethread.h"
ExtractBackground::ExtractBackground(
        std::vector<QWidget*>&inputs,
        QTimer& gl_timer,
        std::vector<arma::uvec>&labels,
        QWidget *parent
        ):
    QFrame(parent),
    inputs_(inputs),
    gl_timer_(gl_timer),
    labels_(labels),
    ui(new Ui::ExtractBackground)
{
    ui->setupUi(this);
    connect(ui->extract_n_plane,SIGNAL(clicked(bool)),this,SLOT(reset_extract_planes()));
    connect(ui->extract_n_plane,SIGNAL(clicked(bool)),this,SLOT(start_extract_planes()));
}

bool ExtractBackground::configure(Config::Ptr)
{
    if(inputs_.empty())return false;
    if(labels_.size()!=inputs_.size())return false;
    return true;
}

ExtractBackground::~ExtractBackground()
{
    delete ui;
}

void ExtractBackground::start_extract_planes()
{
    MeshPairViewerWidget* w = (MeshPairViewerWidget*)(inputs_[current_frame_]);
    ExtractPlaneThread* th = new ExtractPlaneThread(w->first_ptr(),labels_[current_frame_]);
    th->setPlaneNumber(ui->n_planes->value());
    connect(th,SIGNAL(finished()),this,SLOT(finish_extract_planes()));
    th->start(QThread::NormalPriority);
}

void ExtractBackground::finish_extract_planes()
{
    ExtractPlaneThread* th = qobject_cast<ExtractPlaneThread*>(sender());
    if(!th)std::cerr<<"called by unknown object"<<std::endl;
    else{
       th->deleteLater();
    }
    ++current_frame_;
    if(current_frame_<inputs_.size())
    {
        start_extract_planes();
    }
}

void ExtractBackground::extract_points()
{
    ;
}

void ExtractBackground::unextract_points()
{
    ;
}
