#include "regiongrowrgbthread.h"
#include "featurecore.h"
#include "segmentationcore.h"
#include <QTime>
RegionGrowRGBThread::RegionGrowRGBThread(
        std::vector<MeshBundle<DefaultMesh>::Ptr> &inputs,
        std::vector<arma::uvec> &outputs,
        QObject *parent) :
    inputs_(inputs),
    labels_(outputs),
    QObject(parent)
{
    ;
}

bool RegionGrowRGBThread::configure(Config::Ptr config)
{
    config_ = config;
    if(config_->has("JRCS_verbose"))
    {
        verbose_=config_->getInt("JRCS_verbose");
    }
    if(!config_->has("NormalEstimation_k"))return false;
    if(!config_->has("RegionGrow_k"))return false;
    if(!config_->has("RegionGrow_r"))return false;
    if(!config_->has("RegionGrow_cluster_min_num"))return false;
    if(!config_->has("RegionGrow_cluster_max_num"))return false;
    if(!config_->has("RegionGrow_max_norm_angle"))return false;
    return true;
}

void RegionGrowRGBThread::process(void)
{
    QTime timer_;
    if(labels_.size()!=inputs_.size())
    {
        labels_.resize(inputs_.size());
    }
    InputIterator iiter;
    OutputIterator oiter;
    Segmentation::RegionGrowingRGB<DefaultMesh> seg;

    seg.setNumberOfNeighbours(config_->getInt("RegionGrow_k"));
    seg.setMinClusterSize(config_->getInt("RegionGrow_cluster_min_num"));
    seg.setMaxClusterSize(config_->getInt("RegionGrow_cluster_max_num"));
    seg.setSmoothnessThreshold(config_->getFloat("RegionGrow_max_norm_angle")*M_PI/180.0);
    seg.setCurvatureThreshold(std::numeric_limits<float>::max());
    seg.setRadiusOfNeighbours(config_->getFloat("RegionGrow_r"));
    seg.setPointColorThreshold(5.0);
    seg.setRegionColorThreshold(6.0);

    if(verbose_>0)std::cerr<<"RegionGrowRGBThread::process(void)"<<std::endl;

    oiter = labels_.begin();
    timer_.restart();
    for(iiter=inputs_.begin();iiter!=inputs_.end();++iiter)
    {
        MeshBundle<DefaultMesh>& input = **iiter;
        emit message("Region Growing On: "+QString::fromStdString(input.name_),0);
        std::shared_ptr<float> curvature;
        if(input.mesh_.has_vertex_normals())
        {
            DefaultMesh mesh;
            mesh.request_vertex_colors();
            mesh.request_vertex_normals();
            mesh = input.mesh_;
            Feature::computePointNormal(mesh,curvature,0.0,config_->getInt("NormalEstimation_k"));
        }else{
            Feature::computePointNormal(input.mesh_,curvature,0.0,config_->getInt("NormalEstimation_k"));
        }
        seg.setInputMesh(&input.mesh_);
        arma::uvec indices = arma::find(*oiter==0);
        if( indices.size() > config_->getInt("RegionGrow_cluster_min_num"))
        {
            seg.setIndices(indices);
            seg.getCurvatures() = curvature;
            seg.extract(*oiter);
        }
//        std::cerr<<"done region grow"<<std::endl;
        //assigning unknown to cloest
        arma::uvec knownIndices = arma::find( *oiter != 0 );
        arma::uvec unknownIndices = arma::find( *oiter == 0 );
        if((!knownIndices.empty())&&!(unknownIndices.empty()))
        {
            arma::fmat data((float*)input.mesh_.points(),3,input.mesh_.n_vertices(),false,true);
//            std::cerr<<"data.n_cols:"<<data.n_cols<<std::endl;
            arma::fmat knownData = data.cols(knownIndices);
            arma::fmat unknownData = data.cols(unknownIndices);
            ArmaKDTreeInterface<arma::fmat> points(knownData);
            KDTreeSingleIndexAdaptor<
                    L2_Simple_Adaptor<float,ArmaKDTreeInterface<arma::fmat>>,
                    ArmaKDTreeInterface<arma::fmat>,
                    3,arma::uword>
                    kdtree(3,points,KDTreeSingleIndexAdaptorParams(1));
            kdtree.buildIndex();
            arma::uvec rindices(1);
            arma::fvec dists(1);
            for(size_t i = 0 ; i < unknownData.n_cols ; ++i )
            {
                kdtree.knnSearch(unknownData.colptr(i),1,rindices.memptr(),dists.memptr());
                (*oiter)(unknownIndices(i)) = (*oiter)(knownIndices(rindices(0)));
            }
        }
        input.custom_color_.fromlabel(*oiter);
        ++oiter;
    }
    QString msg;
    int ms = timer_.elapsed();
    int s = ms/1000;
    ms -= s*1000;
    int m = s/60;
    s -= m*60;
    int h = m/60;
    m -= h*60;
    msg = msg.sprintf("Time Used of Region Grow:%2u:%2u:%2u.%3u",h,m,s,ms);
    emit message(msg,0);
    emit finished();
}
