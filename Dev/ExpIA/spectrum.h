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
    void composeFunc(
            const MeshBundle<DefaultMesh>::Ptr,
            const arma::uword,
            const arma::uvec&
            );
protected slots:
    void saveCoeff(QString fname=QString());
    void saveLambda(QString fname=QString());
    void saveResidue(QString fname=QString());
    void saveReconstructed(QString fname=QString());
    void saveBases(QString fname=QString());
    void saveAll(void);
    void loadFunc(void);
    void updateSpectrum(void);
    void updateFunction(void);
    void updateCentroid(void);
    void evaluateResEnergy(void);
protected:
    double resEnergy(MeshBundle<DefaultMesh>::Ptr,const arma::uword,const arma::vec&);
private:
    Ui::Spectrum *ui;
    Config::Ptr config_;
    MeshList& inputs_;
    Segmentation::NormalizedCuts<DefaultMesh> ncuts_;
    std::vector<std::shared_ptr<arma::vec>> func_;
    std::vector<std::shared_ptr<arma::vec>> re_func_;
    std::vector<std::shared_ptr<arma::mat>> bases_;
    std::vector<std::shared_ptr<arma::vec>> lambda_;
    std::vector<std::shared_ptr<arma::vec>> coeff_;
    std::vector<std::shared_ptr<arma::vec>> res_;
};

#endif // SPECTRUM_H
