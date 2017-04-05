#include "scenemaker.h"
#include "ui_scenemaker.h"

SceneMaker::SceneMaker(
        std::vector<WidgetPtr>& views,
        MeshBundle<DefaultMesh>::PtrList& inputs,
        QWidget *parent
        ):
    mesh_views_(views),
    inputs_(inputs),
    QWidget(parent),
    ui(new Ui::SceneMaker)
{
    ui->setupUi(this);
    setMinimumSize(320,240);
    lst_view_ = new MeshListViewerWidget(this);
    ui->gridLayout_2->addWidget(lst_view_);
}

bool SceneMaker::configure(Config::Ptr)
{
    lst_view_->query_open_file();
    return true;
}

void SceneMaker::save_to_inputs()
{
    ;
}

SceneMaker::~SceneMaker()
{
    delete ui;
}
