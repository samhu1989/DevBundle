#include "convertmain.h"
#include "common.h"
#include <unistd.h>
int convertMain(int argc, char *argv[])
{
    int ch;
    arma::fmat vv;
    arma::fmat vn;
    arma::Mat<uint8_t> vc;
    std::string out_name;
    while( ( ch = getopt(argc,argv,"hm:i:n:c:o:d") ) != -1 )
    {
        switch(ch)
        {
        case 'i':
            vv.load(std::string(optarg),arma::raw_ascii);
            if(vv.n_rows!=3&&vv.n_cols==3)arma::inplace_trans(vv);
            break;
        case 'n':
            vn.load(std::string(optarg),arma::raw_ascii);
            if(vn.n_rows!=3&&vn.n_cols==3)arma::inplace_trans(vn);
            break;
        case 'c':
            vc.load(std::string(optarg),arma::raw_ascii);
            if(vc.n_rows!=3&&vc.n_cols==3)arma::inplace_trans(vc);
            break;
        case 'o':
            out_name = std::string(optarg);
            break;
        default:
            ;
        }
    }
    DefaultMesh mesh;
    if(vv.n_rows==3)
    {
        for(int idx=0;idx<vv.n_cols;++idx)
        {
            mesh.add_vertex(DefaultMesh::Point(vv(0,idx),vv(1,idx),vv(2,idx)));
        }
    }
    if( vn.n_rows == 3 && vn.n_cols == mesh.n_vertices() )
    {
        mesh.request_vertex_normals();
        arma::fmat normal((float*)mesh.vertex_normals(),3,mesh.n_vertices(),false,true);
        normal = vn;
    }
    if( vc.n_rows == 3 && vn.n_cols == mesh.n_vertices() )
    {
        mesh.request_vertex_colors();
        arma::Mat<uint8_t> color((uint8_t*)mesh.vertex_colors(),3,mesh.n_vertices(),false,true);
        color = vc;
    }
    if(!out_name.empty()&& 0 != mesh.n_vertices() )
    {
        OpenMesh::IO::Options opt;
        opt+=OpenMesh::IO::Options::Binary;
        if(mesh.has_vertex_colors())opt+=OpenMesh::IO::Options::VertexColor;
        if(mesh.has_vertex_normals())opt+=OpenMesh::IO::Options::VertexNormal;
        if(!OpenMesh::IO::write_mesh(mesh,out_name,opt,13)){
            std::cerr<<"can't save to: "<<out_name<<std::endl;
        }
    }
    return 0;
}
