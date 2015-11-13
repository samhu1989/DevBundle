#ifndef REGIONGROWTHREAD_H
#define REGIONGROWTHREAD_H
#include <QThread>
#include <armadillo>
#include "MeshType.h"
class RegionGrowThread:public QThread
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr>::iterator InputIterator;
    typedef std::vector<arma::uvec>::iterator OutputIterator;
    RegionGrowThread(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& outputs
            ):inputs_(inputs),labels_(outputs)
    {
        setObjectName("RegionGrowThread");
    }
signals:
    void message(QString,int);
protected:
    void run();
private:
    std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs_;
    std::vector<arma::uvec>& labels_;
};

#endif // REGIONGROWTHREAD_H
