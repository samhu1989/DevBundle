#include "ncut2d.h"
#include <QMessageBox>
NCut2D::NCut2D(QImage& img,arma::uvec& lbl):
input_img_(img),label_(lbl)
{
    setObjectName("NCut2D");
}

bool NCut2D::configure(Config::Ptr config)
{
    return cut_.configure(config);
}

void NCut2D::process(void)
{
    cut_.cutImage(input_img_,label_);
    std::cerr<<input_img_.width()*input_img_.height()<<std::endl;
//    label_.save("./debug/label/01.arma",arma::raw_ascii);
    emit end();
}
