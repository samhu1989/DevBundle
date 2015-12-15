#ifndef EXTRACTPLANETHREAD_H
#define EXTRACTPLANETHREAD_H
#include <QThread>
#include "common.h"
class ExtractPlaneThread:public QThread
{
public:
    ExtractPlaneThread(
            MeshBundle<DefaultMesh>::Ptr input,
            QObject*parent=0
            ):QThread(parent),input_(input){}
protected:
    void run(void);
private:
    MeshBundle<DefaultMesh>::Ptr input_;
};

#endif // EXTRACTPLANETHREAD_H
