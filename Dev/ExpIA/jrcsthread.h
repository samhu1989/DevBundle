#ifndef JRCSTHREAD_H
#define JRCSTHREAD_H
#include "common.h"
#include <QObject>
#include "objectmodel.h"
#include "jrcsbase.h"
#include <QTime>
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
    bool input_extra(
          const MeshBundle<DefaultMesh>::PtrList& inputs
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
    inline int get_obj_n(void){return jrcs_->get_obj_num();}
    inline int get_k(){return jrcs_->evaluate_k();}
    inline void get_rt(JRCS::JRCSBase::TsLst&rt){jrcs_->get_rt(rt);}
    inline void get_lbl(std::vector<arma::uvec>&lbl){jrcs_->get_label(lbl);}
    inline void get_order(std::vector<arma::uvec>&orders){jrcs_->get_order(orders);}
    inline void set_method(std::shared_ptr<JRCS::JRCSBase> method){jrcs_ = method;}
    inline void set_init_method(std::shared_ptr<JRCS::JRCSInitBase> method){jrcs_->set_init_method(method);}
    inline int get_iter(){if(jrcs_)return jrcs_->get_iter_num();else return -1;}
    inline std::string get_method_name()const{if(jrcs_)return jrcs_->name();else return std::string("unkown method");}
public slots:
    void process(void);
    void get_iter_info(void);

signals:
    void message(QString,int);
    void end(void);
private:
    Config::Ptr config_;
    int verbose_;
    QTime t_;
    std::shared_ptr<JRCS::JRCSBase> jrcs_;
};

#endif // JRCSTHREAD_H
