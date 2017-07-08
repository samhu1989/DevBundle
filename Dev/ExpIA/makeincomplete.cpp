#include "makeincomplete.h"
#include <random>
#include <QFileInfo>
#include <QDir>
MakeIncomplete::MakeIncomplete(
        MeshBundle<DefaultMesh>::PtrList inputs,
        QList<float> ratio,
        QString path,
        QObject *parent
        ):inputs_(inputs),ratio_(ratio),path_(path),QThread(parent)
{
    qSort(ratio_.begin(),ratio_.end());
}

void MakeIncomplete::run(void)
{
    for(MeshBundle<DefaultMesh>::PtrList::iterator iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        saveByRatio(*iter);
    }
}

void MakeIncomplete::saveByRatio(MeshBundle<DefaultMesh>::Ptr in)
{
    DefaultMesh& m = in->mesh_;
    int N = m.n_vertices();
    DefaultMesh om;
    std::vector<DefaultMesh::VertexHandle> sortedV;
    int idx = 0;
    foreach(float r,ratio_)
    {
        QString tmp;
        QString currentpath = path_ + tmp.sprintf("/%02d/",idx);
        QDir dir(currentpath);
        if(!dir.exists())
        {
            dir.mkdir("./");
        }
        int Nr = int(std::round(r*float(N)));
        int N_limit = std::min(Nr,int(sortedV.size()));
        for(int i = om.n_vertices(); i < N_limit ; ++i )
        {
            DefaultMesh::VertexHandle vi = sortedV[i];
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
