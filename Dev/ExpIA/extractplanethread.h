#ifndef EXTRACTPLANETHREAD_H
#define EXTRACTPLANETHREAD_H
#include <QThread>
#include "common.h"
class ExtractPlaneThread:public QThread
{
    Q_OBJECT
public:
    ExtractPlaneThread(
            MeshBundle<DefaultMesh>::Ptr input,
            arma::uvec& label,
            QObject*parent=0
            ):QThread(parent),input_(input),output_(label){}
    inline void setPlaneNumber(uint64_t k){k_=k;}
    inline void setThreshold(float th){threshold_=th;}
protected:
    void run(void);
private:
    MeshBundle<DefaultMesh>::Ptr input_;
    arma::uvec& output_;
    uint64_t k_;
    float threshold_;
};

#endif // EXTRACTPLANETHREAD_H
