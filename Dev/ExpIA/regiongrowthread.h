#ifndef REGIONGROWTHREAD_H
#define REGIONGROWTHREAD_H
#include <QThread>
#include <armadillo>
#include "MeshType.h"
class RegionGrowThread:public QThread
{
public:
    RegionGrowThread(
            std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs,
            std::vector<arma::uvec>& outputs
            ):inputs_(inputs),labels_(outputs){}
protected:
    void run();
private:
    std::vector<MeshBundle<DefaultMesh>::Ptr>& inputs_;
    std::vector<arma::uvec>& labels_;
};

#endif // REGIONGROWTHREAD_H
