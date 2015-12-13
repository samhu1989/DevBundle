#include "extractbackground.h"
#include "ui_extractbackground.h"
#include "MeshPairViewerWidget.h"
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
}

bool ExtractBackground::configure(Config::Ptr)
{
    return true;
}

ExtractBackground::~ExtractBackground()
{
    delete ui;
}

void ExtractBackground::extract_planes()
{
    std::vector<QWidget*>::iterator iter;
    size_t index =0;
    for(iter = inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshPairViewerWidget* w = (MeshPairViewerWidget*)(*iter);
        extract_planes(w->first().mesh_,ui->n_planes->value(),labels_[index]);
    }
}

void ExtractBackground::extract_planes(DefaultMesh&mesh,uint32_t k,arma::uvec&labels)
{
    ;
}

void ExtractBackground::extract_points()
{
    ;
}

void ExtractBackground::unextract_points()
{
    ;
}
