#include "jrcsplatedialog.h"
#include "ui_jrcsplatedialog.h"
#include "jrcsprimitive.h"
JRCSPlateDialog::JRCSPlateDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::JRCSPlateDialog),geo_view_(NULL),plate(NULL)
{
    ui->setupUi(this);
    geo_view_ = new MeshPairViewerWidget();
    ui->horizontalLayout->addWidget(geo_view_);
    geo_view_->first_ptr().reset(new MeshBundle<DefaultMesh>());
    geo_view_->second_ptr().reset(new MeshBundle<DefaultMesh>());
    init_plate();
    init_points();
    geo_view_->set_draw_mode("Plate");
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
        { 0.5, 0.5, 0, 0}
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
    ColorArray::colorfromValue((ColorArray::RGB888*)xc.memptr(),xc.n_cols,value);
}

JRCSPlateDialog::~JRCSPlateDialog()
{
    if(geo_view_)geo_view_->deleteLater();
    if(plate)delete plate;
    delete ui;
}
