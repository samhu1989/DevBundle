#ifndef INPATCHGRAPHCUT_H
#define INPATCHGRAPHCUT_H
#include <QThread>
#include "objectmodel.h"
#include "graphcut.h"
#include "common.h"
#include <typeinfo>
#include <QTime>
class InPatchGraphCut:public QThread
{
    Q_OBJECT
public:
    InPatchGraphCut(
            MeshBundle<DefaultMesh>::PtrList&inputmesh,
            std::vector<ObjModel::Ptr>& inputobj,
            std::vector<arma::uvec>& labels,
            QObject* parent=0
            ):QThread(parent),meshes_(inputmesh),objects_(inputobj),labels_(labels)
    {
        setObjectName("InPatchGraphCut");
        current_frame_ = 0;
    }
public:
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
    void sendMatch(int,MeshBundle<DefaultMesh>::Ptr);
protected:
    void run(void);
    void for_each_frame(void);
    void extract_patch_graph(const arma::uvec&);

    bool prepareDataTerm(Segmentation::GraphCut&);
    std::shared_ptr<double> data_;
    std::shared_ptr<DataCost> current_data_;

    bool prepareSmoothTerm(Segmentation::GraphCut&);
    std::shared_ptr<SmoothnessCost> current_smooth_;
    static MRF::CostVal fnCost(int pix1,int pix2,MRF::Label i,MRF::Label j);
    static arma::sp_mat smooth_;
    bool prepareNeighbors(Segmentation::GraphCut&);
    void applyToFrame(const arma::uvec& gc_label,const arma::uvec& label_indices);
protected:
    uint32_t current_frame_;
    uint32_t current_label_;
    DefaultMesh current_patch_mesh_;
    VoxelGraph<DefaultMesh>::Ptr current_patch_graph_;
    uint32_t label_number_;
    uint32_t pix_number_;
private:
    MeshBundle<DefaultMesh>::PtrList& meshes_;
    std::vector<ObjModel::Ptr>& objects_;
    std::vector<arma::uvec>& labels_;
    Config::Ptr config_;
    QTime timer;
};

#endif // INPATCHGRAPHCUT_H
