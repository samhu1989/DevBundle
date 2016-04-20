#include "mergemain.h"
#include <unistd.h>
#include <string>
#include <QDir>
#include <QFileInfo>
#include "common.h"
#include <fstream>
#include <QDebug>
#include <QStringList>
int merge(const std::string& Tfile,const std::string& spath,DefaultMesh& result)
{
    std::ifstream in;
    in.open(Tfile);
    while(!in.eof())
    {
        DefaultMesh mesh;
        std::string file;
        std::string in_name;
        arma::fmat R(3,3);
        arma::fvec t(3);
        in >> file;
        in >> R(0,0); in >> R(0,1); in >> R(0,2); in >> t(0);
        in >> R(1,0); in >> R(1,1); in >> R(1,2); in >> t(1);
        in >> R(2,0); in >> R(2,1); in >> R(2,2); in >> t(2);
        if(file.front()!='#')return -1;
        file.erase(0,1);
        in_name = spath+"/"+file+".ply";
        OpenMesh::IO::Options opt;
        opt+=OpenMesh::IO::Options::Binary;
        opt+=OpenMesh::IO::Options::VertexColor;
        opt+=OpenMesh::IO::Options::VertexNormal;
        mesh.request_vertex_normals();
        mesh.request_vertex_colors();
        if(!OpenMesh::IO::read_mesh(mesh,in_name,opt,13)){
            std::cerr<<"can't load: "<<in_name<<std::endl;
        }
        arma::fmat v((float*)mesh.points(),3,mesh.n_vertices(),false,true);
        v = R*v;
        v.each_col() += t;
        if(mesh.has_vertex_normals())
        {
            arma::fmat n((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
            n = R*n;
        }
        DefaultMesh::VertexIter v_it;
        DefaultMesh::VertexHandle v_handle;
        for(v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
        {
            v_handle = result.add_vertex(mesh.point(*v_it));
            if(mesh.has_vertex_colors())result.set_color(v_handle,mesh.color(*v_it));
            if(mesh.has_vertex_normals())result.set_normal(v_handle,mesh.normal(*v_it));
        }
    }
    return 0;
}

int mergeMain(int argc, char *argv[])
{
    int ch;
    opterr=0;
    std::string path;
    while( ( ch = getopt(argc,argv,"i:") ) !=-1 )
    {
        switch(ch)
        {
        case 'i':
            path = std::string(optarg);
        default:
            ;
        }
    }
    QDir dir_path_;
    QDir transform_path_,source_path_,target_path_;
    dir_path_.setPath(QString::fromStdString(path));
    transform_path_.setPath(dir_path_.absoluteFilePath("./transform/"));
    if(!transform_path_.exists())return -1;
    source_path_.setPath(dir_path_.absoluteFilePath("./source/"));
    if(!source_path_.exists())return -1;
    target_path_.setPath(dir_path_.absoluteFilePath("./target/"));
    if(!target_path_.exists())return -1;
    QStringList namefilter;
    namefilter<<"*.transform";
    QFileInfoList list = transform_path_.entryInfoList(namefilter,QDir::Files|QDir::NoDotAndDotDot,QDir::NoSort);
    if(list.empty())std::cerr<<"list.size():=0"<<std::endl;
    foreach(QFileInfo info,list)
    {
        DefaultMesh mesh;
        mesh.request_vertex_colors();
        mesh.request_vertex_normals();
        if(
           merge(
            info.absoluteFilePath().toStdString(),
            source_path_.absolutePath().toStdString(),
            mesh
           ) != 0
        )return -1;
        std::string out_name = target_path_.absoluteFilePath( info.baseName()+".ply" ).toStdString();
        if(!out_name.empty()&&0!=mesh.n_vertices())
        {
            OpenMesh::IO::Options opt;
            opt+=OpenMesh::IO::Options::Binary;
            if(mesh.has_vertex_colors())opt+=OpenMesh::IO::Options::VertexColor;
            if(mesh.has_vertex_normals())opt+=OpenMesh::IO::Options::VertexNormal;
            if(!OpenMesh::IO::write_mesh(mesh,out_name,opt,13)){
                std::cerr<<"can't save to: "<<out_name<<std::endl;
            }
        }
    }
    return 0;
}
