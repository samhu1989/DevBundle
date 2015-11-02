#ifndef COMPUTENORMALTHREAD_H
#define COMPUTENORMALTHREAD_H
#include <QThread>
#include <MeshType.h>
class ComputeNormalThread:public QThread
{
public slots:
    ComputeNormalThread(DefaultMesh&mesh,QObject* parent=0):
        QThread(parent),
        mesh_(mesh)
    {
        ;
    }
protected:
    void run();
private:
    DefaultMesh& mesh_;
};

#endif // COMPUTENORMALTHREAD_H
