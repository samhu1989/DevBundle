#include "makeincomplete.h"
#include <random>
#include <QFileInfo>
#include <QDir>
MakeIncomplete::MakeIncomplete(
        MeshBundle<DefaultMesh>::PtrList inputs,
        QList<float> ratio,
        QString path,
        QObject *parent
        ):inputs_(inputs),ratio_(ratio),path_(path),QObject(parent)
{
    qSort(ratio_.begin(),ratio_.end());
}

bool MakeIncomplete::configure(Config::Ptr)
{
    return true;
}

void MakeIncomplete::process(void)
{
    std::cerr<<"a"<<std::endl;
    for(MeshBundle<DefaultMesh>::PtrList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        saveByRatio(*iter);
    }
    emit end();
}

typedef std::pair<DefaultMesh::VertexHandle,float> Pair;

bool less(Pair& p0,Pair& p1)
{
    return p0.second < p1.second;
}

void MakeIncomplete::saveByRatio(MeshBundle<DefaultMesh>::Ptr in)
{
    QString msg;
    DefaultMesh& m = in->mesh_;
    msg = tr("processing:")+QString::fromStdString(in->name_);
    emit message(msg,0);
    int N = m.n_vertices();
    DefaultMesh om;
    if(m.has_vertex_normals())om.request_vertex_normals();
    if(m.has_vertex_colors())om.request_vertex_colors();
    std::vector<std::pair<DefaultMesh::VertexHandle,float>> sortedV;
    DefaultMesh::VertexIter viter;
    std::uniform_int_distribution<int> dist(0,N-1);
    int choosen_i = dist(gen_);
    DefaultMesh::Point choosen_p;
    int idx = 0;
    for(viter=m.vertices_begin();viter!=m.vertices_end();++viter)
    {
        if(idx==choosen_i)
        {
            choosen_p = m.point(*viter);
            break;
        }
        ++idx;
    }
    for(viter=m.vertices_begin();viter!=m.vertices_end();++viter)
    {
        DefaultMesh::Point p = m.point(*viter);
        float x = ( choosen_p[0] - p[0] );
        float y = ( choosen_p[1] - p[1] );
        float z = ( choosen_p[2] - p[2] );
        float d = std::sqrt( x*x + y*y + z*z );
        sortedV.emplace_back(*viter,d);
    }
    std::sort(sortedV.begin(),sortedV.end(),less);
    idx = 0;
    foreach(float r,ratio_)
    {
        QString tmp;
        QString currentpath = path_ + tmp.sprintf("/%03d/",int(std::round(r*100.0)));
        QDir dir(currentpath);
        if(!dir.exists())
        {
            dir.mkdir("./");
        }
        int Nr = int(std::round(r*float(N)));
        int N_limit = std::min(Nr,int(sortedV.size()));
        for(int i = om.n_vertices(); i < N_limit ; ++i )
        {
            DefaultMesh::VertexHandle vi = sortedV[i].first;
            DefaultMesh::VertexHandle vo = om.add_vertex(m.point(vi));
            om.set_normal(vo,m.normal(vi));
            om.set_color(vo,m.color(vi));
        }
        OpenMesh::IO::Options opt;
        opt+=OpenMesh::IO::Options::Binary;
        opt+=OpenMesh::IO::Options::VertexColor;
        opt+=OpenMesh::IO::Options::VertexNormal;
        QString filepath = dir.absoluteFilePath(QString::fromStdString(in->name_+".ply"));
        if(!OpenMesh::IO::write_mesh(om,filepath.toStdString(),opt,13)){
            std::cerr<<"can't save to:"<<filepath.toStdString()<<std::endl;
        }
        idx++;
    }
}
