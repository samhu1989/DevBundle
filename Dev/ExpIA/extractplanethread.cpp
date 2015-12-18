#include "extractplanethread.h"
#include "segmentationcore.h"
#include <memory>
void ExtractPlaneThread::run(void)
{
    arma::uvec& labels = output_;
    arma::fvec axis(3,arma::fill::zeros);
    axis(2) = 1.0;
    float eps(M_PI/90.0);
    std::shared_ptr<Segmentation::SegmentationRANSAC> ransac;
    ransac = std::make_shared<Segmentation::SegmentationRANSAC>(Segmentation::SAC_Model::PERPENDICULLAR_PLANE);
    arma::fmat seg_input((float*)input_->mesh_.points(),3,input_->mesh_.n_vertices(),false,true);

    ransac->setAxis(axis);
    ransac->setEpsAngle(eps);
    ransac->setThreshold(threshold_);

    arma::uvec remained = arma::find(labels!=0);
    if(remained.is_empty())labels.fill(1);
    remained = arma::find(labels!=0);
//    std::cerr<<"remained:"<<remained.size()<<std::endl;

    ransac->input(seg_input,remained);
    ransac->extract(labels);

    ransac = std::make_shared<Segmentation::SegmentationRANSAC>(Segmentation::SAC_Model::PARALLEL_PLANE);
    ransac->setAxis(axis);
    ransac->setEpsAngle(eps);
    ransac->setThreshold(threshold_);

    for( size_t index = 0 ; index < k_ - 1 ; ++index )
    {
        remained = arma::find(labels!=0);
//        std::cerr<<"remained:"<<remained.size()<<std::endl;
        if(remained.is_empty())break;
        ransac->input(seg_input,remained);
//        std::cerr<<"extract:"<<remained.size()<<std::endl;
        ransac->extract(labels);
    }
    input_->custom_color_.fromlabel(labels);
}

