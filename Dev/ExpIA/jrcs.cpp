#include "jrcs.h"
#include "jrcsinitexternal.h"
JRCSWork::MatPtrLst JRCSWork::alpha_ptrlst_;
arma::fvec JRCSWork::obj_prob_;
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

void JRCSWork::Init_Bernolli()
{
    alpha_ptrlst_.resize(inputs_.size());

}

bool JRCSWork::optimize(JRCSView* w)
{
    if( alpha_ptrlst_.empty() || obj_prob_.empty() )return false;
    std::shared_ptr<JRCS::JRCSInitBase> init_ptr(new JRCS::JRCSInitExternal(alpha_ptrlst_,obj_prob_));
    w->set_init_method(init_ptr);
    w->start();
    return true;
}
