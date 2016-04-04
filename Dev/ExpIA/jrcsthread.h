#ifndef JRCSTHREAD_H
#define JRCSTHREAD_H
#include "common.h"
#include <QObject>
#include "objectmodel.h"
#include "jrcscore.h"
class JRCSThread:public QObject
{
    Q_OBJECT
public:
    typedef JRCS::JRCSBase::MatPtrLst MatPtrLst;
    typedef JRCS::JRCSBase::CMatPtrLst CMatPtrLst;
    typedef JRCS::JRCSBase::LCMatPtrLst LCMatPtrLst;
    JRCSThread(QObject* parent=0);
    bool configure(Config::Ptr);
    void input(
          const MatPtrLst& vv,
          const MatPtrLst& vn,
          const CMatPtrLst& vc,
          const LCMatPtrLst& vl
         );
    void resetw(
           const MatPtrLst& wv,
           const MatPtrLst& wn,
           const CMatPtrLst& wc
            );
    void setx();
    inline int get_k(){return jrcs_.evaluate_k();}
public slots:
    void process(void);
signals:
    void message(QString,int);
    void end(void);
protected:
private:
    Config::Ptr config_;
    JRCS::JRCSBase jrcs_;
};

#endif // JRCSTHREAD_H
