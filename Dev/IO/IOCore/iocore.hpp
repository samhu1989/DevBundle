#ifndef IOCORE_HPP
#define IOCORE_HPP
#include "iocore.h"
#include "matio.h"
#include <typeinfo>
namespace MATIO{
template<typename _mat>
void save_to_matlab(const _mat& m, const std::string &file, const std::string &var)
{
    void* data = (void*)m.memptr();
    size_t dims[2] = {m.n_rows,m.n_cols};
    mat_t *mat;
    matvar_t *matvar;
    mat = Mat_CreateVer(file.c_str(),NULL,MAT_FT_MAT73);
    std::string varname;
    if(var.empty())
    {
        char tmp[1024];
        if(!_splitpath_s(file.c_str(),NULL,0,NULL,0,tmp,1024,NULL,0))
        {
            varname = std::string(tmp);
        }
        else varname = "X";
    }else varname = var;
    if(typeid(typename _mat::elem_type) == typeid(double))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_DOUBLE,MAT_T_DOUBLE,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(float))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_SINGLE,MAT_T_SINGLE,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(uint8_t))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT8,MAT_T_UINT8,2,dims,data,0);
    }
    Mat_VarWrite( mat, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(mat);
}
}
#endif // IOCORE_HPP
