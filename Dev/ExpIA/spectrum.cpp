#include "spectrum.h"
#include "ui_spectrum.h"

Spectrum::Spectrum(MeshList &inputs, QWidget *parent) :
    inputs_(inputs),
    QFrame(parent),
    ui(new Ui::Spectrum)
{
    ui->setupUi(this);
}

Spectrum::~Spectrum()
{
    delete ui;
}
