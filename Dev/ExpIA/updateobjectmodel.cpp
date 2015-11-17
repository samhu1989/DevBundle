#include "updateobjectmodel.h"
#include "ui_updateobjectmodel.h"

UpdateObjectModel::UpdateObjectModel(IMeshList &inputs, ILabelList &labels, OModelList &outputs, QWidget *parent) :
    QFrame(parent),
    inputs_(inputs),
    labels_(labels),
    outputs_(outputs),
    ui(new Ui::UpdateObjectModel)
{
    ui->setupUi(this);
}

bool UpdateObjectModel::configure(Config::Ptr config)
{
    config_ = config;
    return true;
}

void UpdateObjectModel::startLater()
{
    ;
}

UpdateObjectModel::~UpdateObjectModel()
{
    delete ui;
}
