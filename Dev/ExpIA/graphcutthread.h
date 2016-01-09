#ifndef GRAPHCUTTHREAD_H
#define GRAPHCUTTHREAD_H
#include <QThread>
#include "objectmodel.h"
#include "graphcut.h"
#include "common.h"
#include <typeinfo>
#include <QTime>
#include "nanoflann.hpp"
class GraphCutThread:public QThread
{
    Q_OBJECT
public:
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<float,MeshKDTreeInterface<DefaultMesh>>,
            MeshKDTreeInterface<DefaultMesh>,
            3,arma::uword> MeshTree;
    typedef nanoflann::KDTreeSingleIndexAdaptor<
            nanoflann::L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
            ArmaKDTreeInterface<arma::fmat>,
            3,arma::uword> ArmaTree;
    typedef MeshKDTreeInterface<DefaultMesh> MTInterface;
    typedef ArmaKDTreeInterface<arma::fmat> ATInterface;
    GraphCutThread(
            MeshBundle<DefaultMesh>::PtrList&inputmesh,
            std::vector<ObjModel::Ptr>& inputobj,
            std::vector<arma::uvec>& outputlabels,
            QObject* parent=0
            ):QThread(parent),meshes_(inputmesh),objects_(inputobj),outputs_(outputlabels)
    {
        setObjectName("GraphCutThread");
        current_frame_ = 0;
    }

public:
    bool configure(Config::Ptr config);
signals:
    void message(QString,int);
    void sendMatch(int,MeshBundle<DefaultMesh>::Ptr);
protected:
    void run(void);
    void showMatch(size_t,DefaultMesh&);
    void showData(size_t);
    void showSmooth();

    /*match object model to get data term*/
    bool prepareDataTermLabelWise();
    void prepareDataForLabel(
            uint32_t l,
            VoxelGraph<DefaultMesh>& graph,
            DefaultMesh& obj,
            arma::fvec &dist_score,
            arma::fvec &norm_score,
            arma::fvec &color_score
            );
    void prepareDataForUnknown();
    void normalizeData();

    /*match transformation through object model to data term to get data term*/
    bool prepareDataTermPixWise();
    void prepareDataForPix(uint32_t,arma::mat&);
    void matchPixtoObject(
            uint32_t pix,
            uint32_t objIdx,
            const arma::fmat &R,
            const arma::fvec &t,
            double& score
            );
    void matchPixtoFrame(
            uint32_t pix,
            uint32_t frameIdx,
            const arma::fmat &R,
            const arma::fvec &t,
            double& score
         );

    std::shared_ptr<double> data_;
    std::shared_ptr<DataCost> current_data_;

    bool prepareSmoothTerm(Segmentation::GraphCut&);
    std::shared_ptr<SmoothnessCost> current_smooth_;
    static MRF::CostVal fnCost(int pix1,int pix2,MRF::Label i,MRF::Label j);
    static arma::sp_mat smooth_;

    bool prepareNeighbors(Segmentation::GraphCut&);

protected:
    uint32_t current_frame_;
    uint32_t label_number_;
    uint32_t pix_number_;
    std::vector<std::shared_ptr<ArmaTree>> mesh_trees_;
    std::vector<std::shared_ptr<ATInterface>> mesh_tree_interface_;
    std::vector<std::shared_ptr<MeshTree>> obj_trees_;
    std::vector<std::shared_ptr<MTInterface>> obj_tree_interface_;
private:
    MeshBundle<DefaultMesh>::PtrList& meshes_;
    std::vector<ObjModel::Ptr>& objects_;
    std::vector<arma::uvec>& outputs_;
    Config::Ptr config_;
    QTime timer;
};

#endif // GRAPHCUTTHREAD_H
