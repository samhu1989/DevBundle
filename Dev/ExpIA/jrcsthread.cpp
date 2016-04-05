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
    }
    return true;
}

void JRCSThread::input(
      const MatPtrLst& vv,
      const MatPtrLst& vn,
      const CMatPtrLst& vc,
      const LCMatPtrLst& vl
     )
{
    jrcs_.input(vv,vn,vc,vl);
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
}

void JRCSThread::process(void)
{
    jrcs_.compute();
    QThread::sleep(8);
    emit end();
}
