#include "sort_agd.h"
#include "agd.h"
Sort_AGD::Sort_AGD(
        MeshBundle<DefaultMesh>::PtrList& inputs,
        QObject *parent
        ):QObject(parent),inputs_(inputs)
{
    ;
}

bool Sort_AGD::configure(Config::Ptr)
{
    if(inputs_.empty()){
        std::cerr<<"Empty Input"<<std::endl;
        return false;
    }
    return true;
}

void Sort_AGD::process()
{
    InputList::iterator iter;
    Feature::AGD<DefaultMesh> agd;
    for(iter=inputs_.begin();iter!=inputs_.end();++iter)
    {
        MeshBundle<DefaultMesh>& m = **iter;
        arma::vec agd_vec;
        agd.extract(m.graph_,agd_vec);
        sort(agd_vec,m);
    }
    emit message(tr("Sort_AGD::process() is finishing"),-1);
    emit finished();
}

void Sort_AGD::sort(const arma::vec& agd,MeshBundle<DefaultMesh>& m)
{
    ;
}
