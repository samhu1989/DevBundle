#include "jrcsthread.h"
#include <QThread>
#include <QTimer>
JRCSThread::JRCSThread(QObject* parent):QObject(parent),jrcs_(new JRCS::JRCSBase())
{
    ;
}

bool JRCSThread::configure(Config::Ptr config)
{
    config_ = config;
    if(config_->has("JRCS_verbose"))
    {
        verbose_=config_->getInt("JRCS_verbose");
    }else{
        verbose_=-1;
    }
    return jrcs_->configure(config);
}

void JRCSThread::input(
      const MatPtrLst& vv,
      const MatPtrLst& vn,
      const CMatPtrLst& vc,
      const LCMatPtrLst& vl
     )
{
    emit message(tr("JRCSThread::input"),0);
    jrcs_->input(vv,vn,vc,vl);
}

bool JRCSThread::input_extra(
      const MeshBundle<DefaultMesh>::PtrList& inputs
     )
{
    return jrcs_->input_extra(inputs);
}

void JRCSThread::input_with_label(
      const MatPtrLst& vv,
      const MatPtrLst& vn,
      const CMatPtrLst& vc,
      const LCMatPtrLst& vlc,
      const LMatPtrLst& vl
     )
{
    emit message(tr("JRCSThread::input_with_label"),0);
    jrcs_->input_with_label(vv,vn,vc,vlc,vl);
}

void JRCSThread::resetw(
       const MatPtrLst& wv,
       const MatPtrLst& wn,
       const CMatPtrLst& wc
        )
{
    jrcs_->resetw(wv,wn,wc);
}

void JRCSThread::resetx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{

    jrcs_->initx(xv,xn,xc);
    jrcs_->reset_rt();
    if(verbose_>0)std::cerr<<"Done resetx"<<std::endl;
}

void JRCSThread::process(void)
{
    t_.restart();
    if(verbose_>0)std::cerr<<"Running "<<jrcs_->name()<<std::endl;
    jrcs_->reset_iteration();
    jrcs_->compute();
    emit end();
}

void JRCSThread::get_iter_info()
{
    QString msg;
    int n = jrcs_->get_max_iter();
    int init_n = jrcs_->get_max_init_iter();
    int iter_i = jrcs_->get_iter_num();
    if(init_n==0)msg = msg.sprintf("iter:%u/%u,in %d ms",iter_i,n,t_.elapsed());
    else {
        if( iter_i <= init_n )
        {
            msg = msg.sprintf("init:%u/%u,in %d ms",iter_i,init_n,t_.elapsed());
        }else{
            msg = msg.sprintf("iter:%u/%u,in %d ms",iter_i,n,t_.elapsed());
        }
    }
    emit message(msg,0);
}
