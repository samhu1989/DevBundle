#include "jrcsthread.h"
#include <QThread>
#include <QTimer>
JRCSThread::JRCSThread(QObject* parent):QObject(parent)
{
    ;
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
    }else return false;

    if(config_->has("JRCS_max_init"))
    {
        jrcs_.set_max_init_iter(config_->getInt("JRCS_max_init"));
    }else jrcs_.set_max_init_iter(config_->getInt("JRCS_max_iter")/2);

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
            if(config_->has("JRCS_smooth_iter"))
            {
                jrcs_.set_max_smooth_iter(config_->getInt("JRCS_smooth_iter"));
            }else jrcs_.set_max_smooth_iter(1);
        }
    }else{
        jrcs_.enable_smooth(true);
        if(config_->has("JRCS_smooth_w"))
        {
            jrcs_.set_smooth_weight(config_->getFloat("JRCS_smooth_w"));
        }else jrcs_.set_smooth_weight(1.0);
        if(config_->has("JRCS_smooth_iter"))
        {
            jrcs_.set_max_smooth_iter(config_->getInt("JRCS_smooth_iter"));
        }else jrcs_.set_max_smooth_iter(1);
    }

    if(config_->has("JRCS_debug_path"))
    {
        jrcs_.set_debug_path(config_->getString("JRCS_debug_path"));
    }else jrcs_.set_debug_path("./debug/");

    if(config_->has("JRCS_mu_type"))
    {
        if(config_->getString("JRCS_mu_type")=="ObjOnly")jrcs_.set_mu_type(JRCS::JRCSBase::ObjOnly);
        if(config_->getString("JRCS_mu_type")=="ObjPointDist")jrcs_.set_mu_type(JRCS::JRCSBase::ObjPointDist);
    }else jrcs_.set_mu_type(JRCS::JRCSBase::ObjOnly);

    if(config_->has("JRCS_rt_type"))
    {
        if(config_->getString("JRCS_rt_type")=="Gamma")jrcs_.set_rt_type(JRCS::JRCSBase::Gamma);
    }else jrcs_.set_rt_type(JRCS::JRCSBase::All);

    return true;
}

void JRCSThread::input(
      const MatPtrLst& vv,
      const MatPtrLst& vn,
      const CMatPtrLst& vc,
      const LCMatPtrLst& vl
     )
{
    emit message(tr("JRCSThread::input"),0);
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

    jrcs_.initx(xv,xn,xc);
    jrcs_.reset_rt();
}

void JRCSThread::process(void)
{
    jrcs_.reset_iteration();
    jrcs_.compute();
    emit end();
}

void JRCSThread::get_iter_info()
{
    QString msg;
    msg = msg.sprintf("iter:%u/%u",jrcs_.get_iter_num(),jrcs_.get_max_iter());
    emit message(msg,0);
}
