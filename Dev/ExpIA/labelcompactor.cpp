#include "labelcompactor.h"
#include <QTime>
LabelCompactor::LabelCompactor(
        MeshList& inputs,
        LabelList& labels,
        QObject *parent):inputs_(inputs),labels_(labels),QObject(parent)
{
    std::cerr<<"LabelCompactor"<<std::endl;
}

bool LabelCompactor::configure(Config::Ptr)
{
    if(labels_.empty())return false;
    if(inputs_.empty())return false;
    if(labels_.size()!=inputs_.size())return false;
    return true;
}

void LabelCompactor::process()
{
    QTime timer_;
    timer_.restart();
    LabelList::iterator iter;
    MeshList::iterator iiter = inputs_.begin();
    QString msg;
    arma::uword f=0;
    for(iter=labels_.begin();iter!=labels_.end();++iter)
    {
        msg = msg.sprintf("Processing Frame %u",f);
        ++f;
        emit message(msg,0);
        arma::uvec& label = *iter;
        arma::uword max = arma::max(label);
        arma::uword cnt = 0;
        for(arma::uword index=0; index <= max ; ++index )
        {
            arma::uvec indices = arma::find(label==index);
            if(!indices.empty()){
                label(indices).fill(cnt);
                ++cnt;
            }
        }
        MeshBundle<DefaultMesh>& m = **iiter;
        m.custom_color_.fromlabel(label);
        ++iiter;
        if(iiter==inputs_.end())break;
    }
    int ms = timer_.elapsed();
    int s = ms/1000;
    ms -= s*1000;
    int m = s/60;
    s -= m*60;
    int h = m/60;
    m -= h*60;
    msg = msg.sprintf("Time Used of Compacting:%02u:%02u:%02u.%03u",h,m,s,ms);
    emit message(msg,0);
    emit end();
}
