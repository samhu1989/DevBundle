#include "jrcsthread.h"
#include <QThread>
JRCSThread::JRCSThread(QObject* parent):QObject(parent)
{

}

bool JRCSThread::configure(Config::Ptr config)
{
    config_ = config;

    if(config_->has("JRCS_obj_w"))
    {
        std::vector<float> objw;
        config_->getFloatVec("JRCS_obj_w",objw);
        jrcs_.reset_objw(objw);
    }else return false;

    if(config_->has("JRCS_max_iter"))
    {
        jrcs_.set_max_iter(config_->getInt("JRCS_max_iter"));
    }
    else return false;

    if(config_->has("JRCS_verbose"))
    {
        if(!config_->getInt("JRCS_verbose") )verbose_=false;
        else verbose_=true;
    }else{
        verbose_=true;
    }

    if(config_->has("JRCS_smooth"))
    {
        if(!config_->getInt("JRCS_smooth"))jrcs_.enable_smooth(false);
        else{
            jrcs_.enable_smooth(true);
            if(config_->has("JRCS_smooth_w"))
            {
                jrcs_.set_smooth_weight(config_->getFloat("JRCS_smooth_w"));
            }else jrcs_.set_smooth_weight(1.0);
        }
    }else{
        jrcs_.enable_smooth(true);
        if(config_->has("JRCS_smooth_w"))
        {
            jrcs_.set_smooth_weight(config_->getFloat("JRCS_smooth_w"));
        }else jrcs_.set_smooth_weight(1.0);
    }

    if(config_->has("JRCS_debug_path"))
    {
        jrcs_.set_debug_path(config_->getString("JRCS_debug_path"));
    }else jrcs_.set_debug_path("./debug/");
    return true;
}

void JRCSThread::input(
      const MatPtrLst& vv,
      const MatPtrLst& vn,
      const CMatPtrLst& vc,
      const LCMatPtrLst& vl
     )
{
    jrcs_.input(vv,vn,vc,vl,verbose_);
}

void JRCSThread::resetw(
       const MatPtrLst& wv,
       const MatPtrLst& wn,
       const CMatPtrLst& wc
        )
{
    jrcs_.resetw(wv,wn,wc);
}

void JRCSThread::resetx(
        const MatPtr& xv,
        const MatPtr& xn,
        const CMatPtr& xc
        )
{
    std::cerr<<"JRCSThread::resetx"<<std::endl;
    jrcs_.initx(xv,xn,xc);
    jrcs_.reset_rt();
}

void JRCSThread::process(void)
{
    jrcs_.reset_iteration();
    jrcs_.compute();
    emit end();
}
