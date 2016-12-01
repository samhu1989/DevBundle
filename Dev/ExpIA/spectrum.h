#ifndef SPECTRUM_H
#define SPECTRUM_H
#include <QFrame>
#include "common.h"
#include "normalizedcuts.h"
namespace Ui {
class Spectrum;
}

class Spectrum : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    explicit Spectrum(
            MeshList& inputs,
            Config::Ptr config,
            QWidget *parent = 0
            );
    ~Spectrum();
protected:
    void updateBase();
    void updateFunc();
    void toPixFunc(MeshBundle<DefaultMesh>::Ptr m,const arma::vec& voxFunc,arma::vec& pixFunc);
    void toVoxFunc(MeshBundle<DefaultMesh>::Ptr m,const arma::vec& pixFunc,arma::vec& voxFunc);
protected slots:
    void saveCoeff(void);
    void saveLambda(void);
    void loadFunc(void);
    void updateSpectrum(void);
    void updateFunction(void);
private:
    Ui::Spectrum *ui;
    Config::Ptr config_;
    MeshList& inputs_;
    Segmentation::NormalizedCuts<DefaultMesh> ncuts_;
    std::vector<std::shared_ptr<arma::vec>> func_;
    std::vector<std::shared_ptr<arma::mat>> bases_;
    std::vector<std::shared_ptr<arma::vec>> lambda_;
    std::vector<std::shared_ptr<arma::vec>> coeff_;
};

#endif // SPECTRUM_H
