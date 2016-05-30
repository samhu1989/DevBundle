#ifndef JRCSTHREAD_H
#define JRCSTHREAD_H
#include "common.h"
#include <QObject>
#include "objectmodel.h"
#include "jrcsbase.h"
class JRCSThread:public QObject
{
    Q_OBJECT
public:
    typedef JRCS::JRCSBase::MatPtr MatPtr;
    typedef JRCS::JRCSBase::CMatPtr CMatPtr;
    typedef JRCS::JRCSBase::MatPtrLst MatPtrLst;
    typedef JRCS::JRCSBase::CMatPtrLst CMatPtrLst;
    typedef JRCS::JRCSBase::LCMatPtrLst LCMatPtrLst;
    typedef JRCS::JRCSBase::LMatPtr LMatPtr;
    typedef JRCS::JRCSBase::LMatPtrLst LMatPtrLst;
    JRCSThread(QObject* parent=0);
    bool configure(Config::Ptr);
    void input(
          const MatPtrLst& vv,
          const MatPtrLst& vn,
          const CMatPtrLst& vc,
          const LCMatPtrLst& vl
         );
    void input_with_label(
          const MatPtrLst& vv,
          const MatPtrLst& vn,
          const CMatPtrLst& vc,
          const LCMatPtrLst& vlc,
          const LMatPtrLst& vl
         );
    void resetw(
           const MatPtrLst& wv,
           const MatPtrLst& wn,
           const CMatPtrLst& wc
            );
    void resetx(
            const MatPtr& xv,
            const MatPtr& xn,
            const CMatPtr& xc
            );
    inline int get_obj_n(void){return jrcs_.get_obj_num();}
    inline int get_k(){return jrcs_.evaluate_k();}
    inline void get_rt(JRCS::JRCSBase::TsLst&rt){jrcs_.get_rt(rt);}
    inline void get_lbl(std::vector<arma::uvec>&lbl){jrcs_.get_label(lbl);}
public slots:
    void process(void);
    void get_iter_info(void);
signals:
    void message(QString,int);
    void end(void);
private:
    Config::Ptr config_;
    int verbose_;
    JRCS::JRCSBase jrcs_;
};

#endif // JRCSTHREAD_H
