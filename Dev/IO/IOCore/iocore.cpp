#include "iocore.h"
#include "iocore.hpp"
#include "iocore_global.h"
#include <armadillo>
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::mat>(const arma::mat&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::vec>(const arma::vec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::rowvec>(const arma::rowvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::fmat>(const arma::fmat&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::fvec>(const arma::fvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::frowvec>(const arma::frowvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT MATIO::save_to_matlab<arma::Mat<uint8_t>>(const arma::Mat<uint8_t>&,const std::string&,const std::string&);
