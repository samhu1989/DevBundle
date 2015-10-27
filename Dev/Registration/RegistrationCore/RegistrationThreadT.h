#ifndef REGISTRATIONTHREAD_H
#define REGISTRATIONTHREAD_H
#include <vector>
#include <QThread>
#include <MeshType.h>
namespace Registration{
    template <typename Reg,typename M>
    class RegistrationThreadT:public QThread
    {
    public:
        typedef typename std::vector<typename MeshBundle<M>::Ptr> MeshList;
        RegistrationThreadT(QObject* parent=0);
        virtual bool init(MeshList& mesh_list);
        const std::string& errorString(void){ return reg_.errorForThread(); }
    protected:
        virtual void compute(void);
    private:
        Reg reg_;
        std::shared_ptr<void*> info_;
    };
}
#include <RegistrationThreadT.hpp>
#endif // REGISTRATIONTHREAD_H
