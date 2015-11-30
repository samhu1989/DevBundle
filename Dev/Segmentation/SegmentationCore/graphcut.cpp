#include "graphcut.h"

#include <strstream>
namespace Segmentation{
GraphCut::GraphCut()
{

}

GraphCut::~GraphCut()
{

}

void GraphCut::setLabelNumber(int n)
{
    numberofLabels = n ;
}

void GraphCut::inputDataTerm(std::shared_ptr<DataCost> data)
{
    data_ = data;
}

void GraphCut::inputSmoothTerm(std::shared_ptr<SmoothnessCost> smooth)
{
    smooth_ = smooth;
}

void GraphCut::optimize(uint32_t iter,float& time)
{
    mrf_->initialize();
    mrf_->clearAnswer();
    mrf_->optimize(iter,time);
    std::stringstream stream;
    stream<<"GraphCut:"
          <<"Done "<<iter<<"iterations in"<<time<<"s"
          <<"Result in Data("<<mrf_->dataEnergy()<<")+"
          <<"Smooth("<<mrf_->smoothnessEnergy()<<")";
    info_ = stream.str();
}

void GraphCut::getAnswer(arma::uvec&label)
{
    label = arma::uvec(numberofPixels);
    MRF::Label* answer = mrf_->getAnswerPtr();
    arma::Col<MRF::Label> l(answer,numberofPixels,false,true);
    label = arma::conv_to<arma::uvec>::from(l);
}

void GraphCut::init(VoxelGraph<DefaultMesh>& graph,Method method)
{
    numberofPixels = graph.voxel_centers.n_cols;
    numberofLabels = 1 + arma::max(graph.voxel_label);
    eng_.reset(new EnergyFunction(data_.get(),smooth_.get()));
    switch(method)
    {
    case EXPANSION:
        mrf_.reset(new Expansion(numberofPixels,numberofLabels,eng_.get()));
        break;
    case  SWAP:
        mrf_.reset(new Swap(numberofPixels,numberofLabels,eng_.get()));
        break;
    case  BELIEF:
        mrf_.reset(new MaxProdBP(numberofPixels,numberofLabels,eng_.get()));
        break;
    }
    std::stringstream stream;
    stream<<"GraphCut:"
          <<"Init with Data("<<mrf_->dataEnergy()<<")+"
          <<"Smooth("<<mrf_->smoothnessEnergy()<<")";
    info_ = stream.str();
}

void GraphCut::setNeighbors(int pix1, int pix2, MRF::CostVal w)
{
    mrf_->setNeighbors(pix1,pix2,w);
}

}
