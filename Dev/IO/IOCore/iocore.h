#ifndef IOCORE_H
#define IOCORE_H
#include <string>
namespace MATIO {
template<typename _mat>
void save_to_matlab(const _mat& m,const std::string& file,const std::string& var=std::string(""));
}
#endif // IOCORE_H
