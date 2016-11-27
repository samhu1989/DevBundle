#ifndef SPECTRUM_H
#define SPECTRUM_H
#include <QFrame>
#include "common.h"
namespace Ui {
class Spectrum;
}

class Spectrum : public QFrame
{
    Q_OBJECT
public:
    typedef std::vector<MeshBundle<DefaultMesh>::Ptr> MeshList;
    explicit Spectrum(MeshList& inputs,QWidget *parent = 0);
    ~Spectrum();
private:
    Ui::Spectrum *ui;
    MeshList& inputs_;
};

#endif // SPECTRUM_H
