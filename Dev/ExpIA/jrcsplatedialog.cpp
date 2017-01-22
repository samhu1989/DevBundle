#include "jrcsplatedialog.h"
#include "ui_jrcsplatedialog.h"
#include "jrcsprimitive.h"
JRCSPlateDialog::JRCSPlateDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::JRCSPlateDialog),geo_view_(NULL),plate(NULL),timer_(new QTimer()),dt_(M_PI/180.0/500.0),dR_(3,3,arma::fill::eye)
{
    ui->setupUi(this);
    geo_view_ = new MeshPairViewerWidget();
    ui->horizontalLayout->addWidget(geo_view_);
    geo_view_->first_ptr().reset(new MeshBundle<DefaultMesh>());
    geo_view_->second_ptr().reset(new MeshBundle<DefaultMesh>());
    init_plate();
    init_points();
    geo_view_->set_draw_mode("Plate");
    geo_view_->set_center_at_mesh(geo_view_->first().mesh_);
    geo_view_->set_normal_scale(0.2);
    connect(ui->tranform,SIGNAL(clicked(bool)),this,SLOT(start_transform()));
    connect(ui->fit,SIGNAL(clicked(bool)),this,SLOT(start_fit()));
    timer_->setSingleShot(false);
    timer_->setInterval(100);
    float dtheta = M_PI / 180 * 1.0;
    arma::fmat Ra = {{std::cos(dtheta),std::sin(dtheta),0},{-std::sin(dtheta),std::cos(dtheta),0},{0,0,1}};
    arma::fmat Rb = {{1,0,0},{0,std::cos(dtheta),std::sin(dtheta)},{0,-std::sin(dtheta),std::cos(dtheta)}};
    arma::fmat Rc = {{std::cos(dtheta),0,-std::sin(dtheta)},{0,1,0},{std::sin(dtheta),0,std::cos(dtheta)}};
    dR_ = Ra*Rb*Rc;
}

void JRCSPlateDialog::init_plate()
{
    DefaultMesh& mesh = geo_view_->first_ptr()->mesh_;

    std::vector<DefaultMesh::VertexHandle>  face_vhandles_a,face_vhandles_b;
    face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
    face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
    face_vhandles_a.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
    face_vhandles_b.push_back(mesh.add_vertex(DefaultMesh::Point(0,0,0)));
    face_vhandles_b.push_back(face_vhandles_a[0]);
    face_vhandles_b.push_back(face_vhandles_a[2]);

    mesh.add_face(face_vhandles_a);
    mesh.add_face(face_vhandles_b);

    mesh.request_face_normals();
    mesh.request_face_colors();
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();

    arma::fmat v = {
        {-0.5, 0.5, 0.5,-0.5},
        { 0, 0, 0, 0},
        { 0.5, 0.5, -0.5, -0.5}
    };
    arma::fmat n = {
        { 0, 0, 0, 0},
        { 1, 1, 1, 1},
        { 0, 0, 0, 0}
    };
    arma::Mat<uint8_t> c = {
        {137,137,137,137},
        {157,157,157,157},
        {192,192,192,192}
    };

    arma::fmat xv((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::fmat xn((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
    arma::Mat<uint8_t> xc((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);

    xv = v; xn = n; xc = c;

    arma::fvec pos = {0,0,0};

    plate = new JRCS::Plate(xv,xn,xc,pos);
}

void JRCSPlateDialog::init_points()
{
    DefaultMesh& mesh = geo_view_->second_ptr()->mesh_;
    int pN = 1000;
    int i = 0;
    while( i < pN )
    {
        mesh.add_vertex(DefaultMesh::Point(0,0,0));
        ++i;
    }
    mesh.request_vertex_normals();
    mesh.request_vertex_colors();
    arma::fmat xv((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::fmat xn((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
    arma::Mat<uint8_t> xc((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);

    arma::vec value(mesh.n_vertices());
    i = 0;
    while( i < pN )
    {
        arma::fmat tmp(3,100,arma::fill::randu);
        tmp -= 0.5;
        tmp *= 2.0;
        arma::vec dist = plate->get_dist2(tmp);
        arma::uword idx = 0;
        value(i) = dist.min(idx);
        xv.col(i) = tmp.col(idx);
        ++ i;
    }
    xn.fill(0.0);
    xn.row(1).fill(1.0);
    ColorArray::colorfromValue((ColorArray::RGB888*)xc.memptr(),xc.n_cols,arma::sqrt(value));
}

arma::fvec helix(float t)
{
    arma::fvec pos = arma::fvec(3,arma::fill::zeros);
    pos(0) = 1.0*std::sin(180.0*t)*std::cos(360.0*10.0*t);
    pos(1) = 1.0*std::sin(180.0*t)*std::sin(360.0*10.0*t);
    pos(2) = 1.0*std::cos(180.0*t);
    return pos;
}

void JRCSPlateDialog::start_transform()
{
    if(!timer_->isActive());
    {
        time_ = 0.0;
        arma::fvec t = helix( time_ );
        plate->translate(t,*plate);
        translate_ = t;
        connect(timer_,SIGNAL(timeout()),this,SLOT(transform()));
        timer_->start();
    }
}

void JRCSPlateDialog::transform()
{
    timer_->stop();
    DefaultMesh& mesh = geo_view_->second_ptr()->mesh_;
    arma::fmat xv((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::Mat<uint8_t> xc((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);
    arma::vec value(mesh.n_vertices());
    arma::fmat R(3,3,arma::fill::eye);
    arma::fvec t = -translate_;
    plate->translate(t,*plate);
    t = helix(time_);
    R = dR_;
    plate->transform(R,t,*plate);
    translate_ = t;
    arma::fvec scale = {std::pow(0.5,0.002),1.0,std::pow(2.0,0.002)};
    plate->scale(scale,*plate);
    value = plate->get_dist2(xv);
    ColorArray::colorfromValue((ColorArray::RGB888*)xc.memptr(),xc.n_cols,arma::sqrt(value));
    time_ += dt_;
    geo_view_->updateGL();
    if(time_>=500*dt_){
        timer_->disconnect(this,SLOT(transform()));
    }
    else{
        timer_->start();
    }
}

void JRCSPlateDialog::start_fit()
{
    if(!timer_->isActive());
    {
        time_ = 0.0;
        arma::fvec t = {0,0.2,0};
        arma::fvec s = {0.5,0,1.5};
        plate->translate(t,*plate);
        plate->scale(s,*plate);
        translate_ = t;
        connect(timer_,SIGNAL(timeout()),this,SLOT(fit()));
        timer_->start();
    }
}

void JRCSPlateDialog::fit()
{
    timer_->stop();
    std::cerr<<"fitting voked"<<std::endl;
    DefaultMesh& mesh = geo_view_->second_ptr()->mesh_;
    arma::fmat xv((float*)mesh.points(),3,mesh.n_vertices(),false,true);
    arma::fmat xn((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
    arma::Mat<uint8_t> xc((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);
    arma::vec value(mesh.n_vertices());
    value.fill(1.0);
    std::cerr<<"start accumulate"<<std::endl;
    plate->accumulate(
                xv,
                xn,
                xc,
                value
                );
    std::cerr<<"fit"<<std::endl;
    plate->fit();
    ColorArray::colorfromValue((ColorArray::RGB888*)xc.memptr(),xc.n_cols,arma::sqrt(value));
    time_ += dt_;
    geo_view_->updateGL();
    if(time_>=500*dt_){
        timer_->disconnect(this,SLOT(transform()));
    }else{
        timer_->start();
    }
}

JRCSPlateDialog::~JRCSPlateDialog()
{
    if(timer_)delete timer_;
    if(geo_view_)geo_view_->deleteLater();
    if(plate)delete plate;
    delete ui;
}
