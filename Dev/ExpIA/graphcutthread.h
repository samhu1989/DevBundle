#ifndef GRAPHCUTTHREAD_H
#define GRAPHCUTTHREAD_H
#include <QThread>
#include "objectmodel.h"
#include "graphcut.h"
class GraphCutThread:public QThread
{
    Q_OBJECT
public:
    GraphCutThread(
            MeshBundle<DefaultMesh>::PtrList&inputmesh,
            std::vector<ObjModel::Ptr>& inputobj,
            std::vector<arma::uvec>& outputlabels
            ):meshes_(inputmesh),objects_(inputobj),outputs_(outputlabels)
    {
        setObjectName("GraphCutThread");
        current_frame_ = 0;
    }
signals:
    void message(QString,int);
public:
    bool configure(Config::Ptr config);
protected:
    void run(void);

    bool prepareDataTerm();
    bool prepareSmoothTerm();

    uint32_t current_frame_;
    std::shared_ptr<DataCost> current_data_;
    std::shared_ptr<SmoothnessCost> current_smooth_;
private:
    MeshBundle<DefaultMesh>::PtrList& meshes_;
    std::vector<ObjModel::Ptr>& objects_;
    std::vector<arma::uvec>& outputs_;
    Config::Ptr config_;
};

#endif // GRAPHCUTTHREAD_H
