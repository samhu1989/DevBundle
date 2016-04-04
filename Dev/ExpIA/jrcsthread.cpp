#include "jrcsthread.h"
#include <QThread>
JRCSThread::JRCSThread(QObject* parent):QObject(parent)
{

}

bool JRCSThread::configure(Config::Ptr)
{
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

void JRCSThread::process(void)
{
    jrcs_.compute();
    QThread::sleep(2);
}
