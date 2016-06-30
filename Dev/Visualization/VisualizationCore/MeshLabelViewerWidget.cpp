#include "MeshLabelViewerWidget.h"
#include "ui_meshlabelviewerwidget.h"

MeshLabelViewerWidget::MeshLabelViewerWidget(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::MeshLabelViewerWidget)
{
    ui->setupUi(this);
}

MeshLabelViewerWidget::~MeshLabelViewerWidget()
{
    delete ui;
}
