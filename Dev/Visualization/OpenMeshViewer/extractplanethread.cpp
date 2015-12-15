#include "extractplanethread.h"
#include "segmentationcore.h"
#include <memory>
void ExtractPlaneThread::run(void)
{
    arma::fvec axis(3,arma::fill::zeros);
    axis(2) = 1.0;
    float eps(M_PI/90.0);
    std::shared_ptr<Segmentation::SegmentationRANSAC> ransac;
    ransac = std::make_shared<Segmentation::SegmentationRANSAC>(Segmentation::SAC_Model::PERPENDICULLAR_PLANE);
    arma::fmat seg_input((float*)input_->mesh_.points(),3,input_->mesh_.n_vertices(),false,true);
    ransac->setAxis(axis);
    ransac->setEpsAngle(eps);
    ransac->input(seg_input);
    arma::uvec label;
    ransac->extract(label);
    ransac = std::make_shared<Segmentation::SegmentationRANSAC>(Segmentation::SAC_Model::PARALLEL_PLANE);
    ransac->setAxis(axis);
    ransac->setEpsAngle(eps);
    arma::uvec remained = arma::find(label!=0);
    ransac->input(seg_input,remained);
    ransac->extract(label);
    remained = arma::find(label!=0);
    ransac->input(seg_input,remained);
    ransac->extract(label);
    input_->custom_color_.fromlabel(label);
}

