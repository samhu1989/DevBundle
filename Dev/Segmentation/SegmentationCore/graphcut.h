#ifndef GRAPHCUT_H
#define GRAPHCUT_H
#include "segmentationcore_global.h"
#include <armadillo>
#include "common.h"
#include "MRF/include/mrf.h"
#include "MRF/include/GCoptimization.h"
#include "MRF/include/MaxProdBP.h"
namespace Segmentation{
class SEGMENTATIONCORESHARED_EXPORT GraphCut
{
public:
    typedef enum{
        ICM,
        EXPANSION,
        SWAP,
        BELIEF
    }Method;
    GraphCut();
    ~GraphCut();
    void setLabelNumber(size_t n){numberofLabels = n;}
    void setPixelNumber(size_t n){numberofPixels = n;}
    void inputDataTerm(std::shared_ptr<DataCost>);
    void inputSmoothTerm(std::shared_ptr<SmoothnessCost>);
    void init(Method);
    void updateInfo(void);
    bool setNeighbors(int pix1, int pix2, MRF::CostVal w);
    void optimize(uint32_t,float&);
    void getAnswer(arma::uvec&);
    const std::string& info(){return info_;}
private:
    int numberofPixels;
    int numberofLabels;
    std::shared_ptr<DataCost> data_;
    std::shared_ptr<SmoothnessCost> smooth_;
    std::shared_ptr<EnergyFunction> eng_;
    std::shared_ptr<MRF> mrf_;
    std::string info_;
};
}

#endif // GRAPHCUT_H
