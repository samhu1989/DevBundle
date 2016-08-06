#include "ncut2d.h"
NCut2D::NCut2D(QImage& img,arma::uvec& lbl):
input_img_(img),label_(lbl)
{
    setObjectName("NCut2D");
}

void NCut2D::process(void)
{
    arma::sp_mat A = arma::sprandu<arma::sp_mat>(100, 100, 0.1);
    arma::sp_mat B = A.t()*A;
    arma::vec eigval;
    arma::mat eigvec;
    eigs_sym(eigval, eigvec, B, 5);
    cut_.cutImage(input_img_,label_);
    std::cerr<<input_img_.width()*input_img_.height()<<std::endl;
    label_.save("./debug/label/01.arma",arma::raw_ascii);
    emit end();
}
