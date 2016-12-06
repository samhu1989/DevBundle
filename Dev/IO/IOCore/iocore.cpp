#include "iocore.h"
#include "iocore.hpp"
namespace MATIO{
template void IOCORESHARED_EXPORT save_to_matlab<arma::mat>(const arma::mat&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::vec>(const arma::vec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::rowvec>(const arma::rowvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::fmat>(const arma::fmat&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::fvec>(const arma::fvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::frowvec>(const arma::frowvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::Mat<uint8_t>>(const arma::Mat<uint8_t>&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::uvec>(const arma::uvec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::ivec>(const arma::ivec&,const std::string&,const std::string&);
template void IOCORESHARED_EXPORT save_to_matlab<arma::mat>(const std::vector<std::shared_ptr<arma::mat>>&,const std::string&,const std::vector<std::string>&);
template bool IOCORESHARED_EXPORT load_to_arma<arma::uvec>(arma::uvec&,const std::string&,const std::string&);
template bool IOCORESHARED_EXPORT load_to_arma<arma::vec>(arma::vec&,const std::string&,const std::string&);
template bool IOCORESHARED_EXPORT load_to_arma<arma::mat>(arma::mat&,const std::string&,const std::string&);
}
