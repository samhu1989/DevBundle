#ifndef IOCORE_H
#define IOCORE_H
#include "iocore_global.h"
#include <string>
#include <armadillo>
#include <memory>
namespace MATIO {
template<typename _mat>
void save_to_matlab(const _mat& m,const std::string& file,const std::string& var=std::string(""));
template<typename _mat>
void save_to_matlab(const std::vector<std::shared_ptr<_mat>>& m_lst, const std::string &file, const std::vector<std::string> &var);
template<typename _mat>
bool load_to_arma(_mat& m,const std::string& file,const std::string& var=std::string(""));
}
#endif // IOCORE_H
