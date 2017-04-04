#include "scenemaker.h"
#include "ui_scenemaker.h"

SceneMaker::SceneMaker(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::SceneMaker)
{
    ui->setupUi(this);
    lst_view_ = new MeshListViewerWidget(this);
}

bool SceneMaker::configure(Config::Ptr)
{
    return false;
}

SceneMaker::~SceneMaker()
{
    delete ui;
}
