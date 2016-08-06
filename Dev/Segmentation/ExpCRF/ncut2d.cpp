#include "ncut2d.h"

NCut2D::NCut2D(QImage& img,arma::uvec& lbl):
input_img_(img),label_(lbl)
{
    setObjectName("NCut2D");
}

void NCut2D::process(void)
{
    emit end();
}
