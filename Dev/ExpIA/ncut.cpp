#include "ncut.h"
#include <QMessageBox>
NCut::NCut(
        MeshBundle<DefaultMesh>::PtrList& inputs,
        std::vector<arma::uvec> &labels,
        QObject *parent):
    inputs_(inputs),
    labels_(labels),
    QObject(parent)
{
    setObjectName("NCut");
}

bool NCut::configure(Config::Ptr config)
{
    if(inputs_.empty())return false;
    if(inputs_.front()->graph_.empty()){
        QMessageBox::warning(NULL,tr("NCut"),tr("Please build graph by supervoxel first!"));
        return false;
    }
    return cuts_.configure(config);
}

void NCut::process()
{
    std::vector<arma::uvec>::iterator oiter;
    oiter = labels_.begin();
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {

        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        arma::uvec label;
        cuts_.cutGraph(mesh_ptr,label);
        mesh_ptr->graph_.sv2pix(label,*oiter);
        mesh_ptr->custom_color_.fromlabel(*oiter);
        ++oiter;
        if(oiter==labels_.end())break;
    }
    emit end();
}

void NCut::debug_convexity()
{
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        cuts_.debug_convexity(mesh_ptr);
    }
    emit end();
}

void NCut::debug_W()
{
    MeshBundle<DefaultMesh>::PtrList::iterator iiter;
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>::Ptr mesh_ptr = *iiter;
        cuts_.debug_W(mesh_ptr);
    }
    emit end();
}
