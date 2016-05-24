#include "jrcsinitthread.h"
#include "hierarchicalization.h"
JRCSInitThread::JRCSInitThread(
        MeshBundle<DefaultMesh>::PtrList& inputs,
        std::vector<arma::uvec>& labels,
        QObject* parent
        ):inputs_(inputs),labels_(labels),QObject(parent)
{
    ;
}

bool JRCSInitThread::configure(Config::Ptr config)
{
    return seg_.configure(config);
}

void JRCSInitThread::process()
{
    MeshBundle<DefaultMesh>::PtrList::iterator iter;
    uint32_t index = 0;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        QString msg;
        msg = msg.sprintf("JRCSInitThread::process(%u)",index);
        emit message(msg,0);
        MeshBundle<DefaultMesh>& mesh = **iter;
        seg_.compute(mesh.mesh_);
//        std::cerr<<"getting label"<<std::endl;
        seg_.getBBoxLabel(labels_[index]);
//        seg_.getPlaneLabel(labels_[index]);
        MeshBundle<DefaultMesh>::Ptr boxptr(new MeshBundle<DefaultMesh>);
        seg_.getObjectBox(boxptr->mesh_);
        emit showbox(index,boxptr);
        mesh.custom_color_.fromlabel(labels_[index]);
        ++ index;
    }
    emit finished();
}
