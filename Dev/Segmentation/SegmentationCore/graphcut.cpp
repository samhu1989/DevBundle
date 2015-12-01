#include "graphcut.h"

#include <strstream>
namespace Segmentation{
GraphCut::GraphCut()
{

}

GraphCut::~GraphCut()
{

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

void GraphCut::init(Method method)
{
    if(!data_||0==data_.use_count())std::logic_error("!data_||0==data_.use_count()");
    if(!smooth_||0==smooth_.use_count())std::logic_error("!smooth_||0==smooth_.use_count()");
    eng_.reset(new EnergyFunction(data_.get(),smooth_.get()));
    switch(method)
    {
    case EXPANSION:
        std::cerr<<"EXPANSION"<<std::endl;
        mrf_.reset(new Expansion(numberofPixels,numberofLabels,eng_.get()));
        break;
    case  SWAP:
        std::cerr<<"SWAP"<<std::endl;
        mrf_.reset(new Swap(numberofPixels,numberofLabels,eng_.get()));
        break;
    case  BELIEF:
        std::cerr<<"BELIEF"<<std::endl;
        mrf_.reset(new MaxProdBP(numberofPixels,numberofLabels,eng_.get()));
        break;
    default:
        std::cerr<<"EXPANSION"<<std::endl;
        mrf_.reset(new Expansion(numberofPixels,numberofLabels,eng_.get()));
    }
    mrf_->initialize();
    mrf_->clearAnswer();
}

void GraphCut::updateInfo(void)
{
    std::stringstream stream;
    stream<<"GraphCut:";
    std::cerr<<"Get Data Energy"<<std::endl;
    stream<<"Data("<<mrf_->dataEnergy()<<")+";
    std::cerr<<"Get Smoothness Energy"<<std::endl;
    stream<<"Smooth("<<mrf_->smoothnessEnergy()<<")";
    info_ = stream.str();
}

bool GraphCut::setNeighbors(int pix1, int pix2, MRF::CostVal w)
{
    if(!mrf_||0==mrf_.use_count())
    {
        return false;
    }
    else{
        mrf_->setNeighbors(pix1,pix2,w);
        return true;
    }
}

}
