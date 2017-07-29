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
    connect(ui->toolButton_2,SIGNAL(clicked(bool)),this,SLOT(save_rt()));
    connect(lst_view_,SIGNAL(have_been_transfomed(arma::fmat,arma::fvec,int)),this,SLOT(update_rt(arma::fmat,arma::fvec,int)));
}

bool SceneMaker::configure(Config::Ptr)
{
    lst_view_->query_open_file();
    R_lst_.resize(lst_view_->list().size());
    t_lst_.resize(lst_view_->list().size());
    reset_rt();
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
    save_rt();
    reset_rt();
}

void SceneMaker::save_to_input(DefaultMesh& mesh)
{
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter = lst_view_->list().begin();
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator list_end = lst_view_->list().end();
    int idx = -1;
    for(;iter!=list_end;++iter)
    {
        ++idx;
        if( !visible(idx) )continue;
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
    idx = -1;
    for(;iter!=list_end;++iter)
    {
        ++idx;
        if( ! visible(idx) )continue;
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

bool SceneMaker::visible(const int&idx)
{
    uint32_t start = lst_view_->visible_index();
    uint32_t num = lst_view_->visible_num();
    uint32_t total_num = lst_view_->list().size();
    if( idx >= start && idx < std::min(total_num,start+num)) return true;
    if( (start+num) > total_num  && idx >=0 && idx < (start+num)%total_num) return true;
    return false;
}

void SceneMaker::save_to_label(arma::uvec& lbl)
{
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator iter = lst_view_->list().begin();
    std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator list_end = lst_view_->list().end();
    arma::uword startl = 0;
    int idx = -1;
    for(;iter!=list_end;++iter)
    {
        ++idx;
        if( !visible(idx) )continue;
        MeshBundle<DefaultMesh>& bundle = **iter;
        lbl.rows(startl,startl+bundle.mesh_.n_vertices()-1).fill(idx+1);
        startl += bundle.mesh_.n_vertices();
    }
}

void SceneMaker::save_rt()
{
    QString fileName = QFileDialog::getSaveFileName(this,
                                            tr("Save RT"),
                                            tr("./"),
                                            tr("TXT(*.txt)"));
    if(fileName.isEmpty())return;
    std::ofstream out;
    out.open(fileName.toStdString());
    std::vector<arma::fmat>::iterator riter;
    std::vector<arma::fvec>::iterator titer = t_lst_.begin();
    for(riter=R_lst_.begin();riter!=R_lst_.end();++riter)
    {
        arma::fmat& R = *riter;
        arma::fvec& t = *titer;
        out << R(0,0)<<" "<< R(0,1)<<" "<<R(0,2)<<" "<<t(0)<<std::endl;
        out << R(1,0)<<" "<< R(1,1)<<" "<<R(1,2)<<" "<<t(1)<<std::endl;
        out << R(2,0)<<" "<< R(2,1)<<" "<<R(2,2)<<" "<<t(2)<<std::endl;
        ++titer;
        if(titer==t_lst_.end())break;
    }
    out.close();
}

void SceneMaker::reset_rt()
{
    std::vector<arma::fmat>::iterator riter;
    std::vector<arma::fvec>::iterator titer;
    for(riter=R_lst_.begin();riter!=R_lst_.end();++riter)
    {
        *riter = arma::fmat(3,3,arma::fill::eye);
    }
    for(titer=t_lst_.begin();titer!=t_lst_.end();++titer)
    {
        *titer = arma::fvec(3,arma::fill::zeros);
    }
}

void SceneMaker::update_rt(arma::fmat R,arma::fvec t,int idx)
{
    if(R_lst_[idx].empty())
    {
        R_lst_[idx] = R;
    }else
    {
        R_lst_[idx] = R*R_lst_[idx];
    }
    if(t_lst_[idx].empty())
    {
        t_lst_[idx] = t;
    }else
    {
        t_lst_[idx] += ( R*t_lst_[idx] + t );
    }
}

SceneMaker::~SceneMaker()
{
    delete ui;
}
