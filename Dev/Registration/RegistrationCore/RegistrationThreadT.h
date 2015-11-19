#ifndef REGISTRATIONTHREAD_H
#define REGISTRATIONTHREAD_H
#include <vector>
#include <QThread>
#include "common.h"
namespace Registration{
    template <typename Reg,typename M>
    class RegistrationThreadT:public QThread
    {
    public:
        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        RegistrationThreadT(QObject* parent=0);
        virtual bool init(MeshList&,Config::Ptr&);
        virtual bool init(MeshList&,std::vector<arma::uword>&,Config::Ptr&);
        typename Reg::ResPtr result(void);
        const std::string& errorString(void){ return reg_.errorForThread(); }
    protected:
        virtual void compute(void);
    protected:
        Reg reg_;
    };
}
#include <RegistrationThreadT.hpp>
#endif // REGISTRATIONTHREAD_H
