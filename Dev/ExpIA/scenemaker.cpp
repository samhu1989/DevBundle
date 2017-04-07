#include "scenemaker.h"
#include "ui_scenemaker.h"
#include "meshpairviewerwidget.h"
SceneMaker::SceneMaker(
        std::vector<WidgetPtr>& views,
        MeshBundle<DefaultMesh>::PtrList& inputs,
        std::vector<arma::uvec> &labels,
        QWidget *parent
        ):
    mesh_views_(views),
    inputs_(inputs),
    labels_(labels),
    QWidget(parent),
    ui(new Ui::SceneMaker)
{
    ui->setupUi(this);
    setMinimumSize(320,240);
    lst_view_ = new MeshListViewerWidget(this);
    OpenMesh::IO::Options opt;
    opt += OpenMesh::IO::Options::Binary;
    opt += OpenMesh::IO::Options::VertexColor;
    opt += OpenMesh::IO::Options::VertexNormal;
    opt += OpenMesh::IO::Options::VertexTexCoord;
    opt += OpenMesh::IO::Options::FaceColor;
    opt += OpenMesh::IO::Options::FaceNormal;
    opt += OpenMesh::IO::Options::FaceTexCoord;
    lst_view_->setOptions(opt);
    ui->verticalLayout->addWidget(lst_view_);
    connect(ui->toolButton,SIGNAL(clicked(bool)),this,SLOT(save_to_inputs()));
}

bool SceneMaker::configure(Config::Ptr)
{
    lst_view_->query_open_file();
    return true;
}

void SceneMaker::save_to_inputs()
{
    inputs_.push_back(std::make_shared<MeshBundle<DefaultMesh>>());
    QString name;
    name = name.sprintf("%02d",inputs_.size());
    save_to_input(inputs_.back()->mesh_);
    inputs_.back()->name_ = name.toStdString();
    labels_.emplace_back(inputs_.back()->mesh_.n_vertices(),arma::fill::zeros);
    save_to_label(labels_.back());
    inputs_.back()->custom_color_.fromlabel(labels_.back());
    MeshBundle<DefaultMesh>::Ptr& bundle_ptr = inputs_.back();
    MeshPairViewerWidget* widget = new MeshPairViewerWidget(this);
    widget->setMinimumSize(300,200);
    widget->first_ptr() = bundle_ptr;
    widget->set_center_at_mesh(bundle_ptr->mesh_);
    widget->setWindowTitle(QString::fromStdString(bundle_ptr->name_));
    emit view_input(widget);
}

void SceneMaker::save_to_input(DefaultMesh& mesh)
{
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter = lst_view_->list().begin();
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator list_end = lst_view_->list().end();
    for(;iter!=list_end;++iter)
    {
        MeshBundle<DefaultMesh>& bundle = **iter;
        DefaultMesh::VertexIter viter;
        for(viter=bundle.mesh_.vertices_begin();viter!=bundle.mesh_.vertices_end();++viter)
        {
            DefaultMesh::VertexHandle vh = mesh.add_vertex(bundle.mesh_.point(*viter));
        }
    }
    mesh.request_vertex_colors();
    mesh.request_vertex_normals();
    DefaultMesh::VertexIter vviter = mesh.vertices_begin();
    iter = lst_view_->list().begin();
    for(;iter!=list_end;++iter)
    {
        MeshBundle<DefaultMesh>& bundle = **iter;
        DefaultMesh::VertexIter viter;
        for(viter=bundle.mesh_.vertices_begin();viter!=bundle.mesh_.vertices_end();++viter)
        {

            mesh.set_color( *vviter , bundle.mesh_.color(*viter) );
            OpenMesh::Vec3f n = bundle.mesh_.normal(*viter);
            mesh.set_normal(*vviter , n.normalize() );
            ++vviter;
        }
    }
}

void SceneMaker::save_to_label(arma::uvec& lbl)
{
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter = lst_view_->list().begin();
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator list_end = lst_view_->list().end();
    arma::uword start = 0;
    arma::uword oidx = 1;
    for(;iter!=list_end;++iter)
    {
        MeshBundle<DefaultMesh>& bundle = **iter;
        lbl.rows(start,start+bundle.mesh_.n_vertices()-1).fill(oidx);
        ++oidx;
        start += bundle.mesh_.n_vertices();
    }
}

SceneMaker::~SceneMaker()
{
    delete ui;
}
