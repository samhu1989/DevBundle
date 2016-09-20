#include "jrcs.h"
JRCSWork::MatPtrLst JRCSWork::alpha_ptrlst_;
JRCSWork::JRCSWork(
        MeshList& inputs,
        LabelList& labels,
        ModelList& objects,
        QObject *parent
        ) :
    inputs_(inputs),labels_(labels),objects_(objects),
    QObject(parent)
{
    ;
}

bool JRCSWork::configure(Config::Ptr config)
{
    return true;
}

void JRCSWork::Init_SI_HSK()
{
    ;
}
