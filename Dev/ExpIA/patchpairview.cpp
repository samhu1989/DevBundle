#include "patchpairview.h"
#include "ui_patchpairview.h"

PatchPairView::PatchPairView(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::PatchPairView),
    object_view_(new MeshViewerWidget(this)),
    patch_view_(new MeshViewerWidget(this))
{
    ui->setupUi(this);
    ui->layout->addWidget(object_view_);
    ui->layout->addWidget(patch_view_);
}

PatchPairView::~PatchPairView()
{
    delete ui;
}
