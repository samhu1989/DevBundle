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
    matvar_t *matvar = NULL;
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
    if(typeid(typename _mat::elem_type) == typeid(long long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_INT64,MAT_T_INT64,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(unsigned long long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT64,MAT_T_UINT64,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(unsigned long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT32,MAT_T_UINT32,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_INT32,MAT_T_INT32,2,dims,data,0);
    }
    if(!matvar){
        Mat_Close(mat);
        std::cerr<<"MATIO::save_to_matlab(unknown type)"<<std::endl;
        return;
    }
    Mat_VarWrite( mat, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(mat);
}
template<typename _mat>
bool load_to_arma(const _mat& m,const std::string& file,const std::string& var)
{
    mat_t *mat = Mat_Open(file.c_str(),MAT_ACC_RDONLY);
    if(!mat){
        std::cerr<<"MATIO::load_to_arma(failed to load "<<file<<")"<<std::endl;
        return false;
    }
    matvar_t *matvar = NULL;
    std::string varname;
    if(var.empty())
    {
        char tmp[1024];
        if(!_splitpath_s(file.c_str(),NULL,0,NULL,0,tmp,1024,NULL,0))
        {
            varname = std::string(tmp);
        }
        else varname = "X";
    }else {
        varname = var;
    }
    matvar = Mat_VarReadInfo(mat,var.c_str());
    if(!matvar)
    {
        std::cerr<<"MATIO::load_to_arma(failed to find variable "<<varname<<")"<<std::endl;
        return false;
    }
    if(typeid(typename _mat::elem_type) == typeid(double))
    {
        if(matvar->rank!=2)
        {
            std::cerr<<"MATIO::load_to_arma(not a matrix,a data with rank="<<matvar->rank<<")"<<std::endl;
            return false;
        }
        if(matvar->class_type!=MAT_C_DOUBLE||matvar->data_type!=MAT_T_DOUBLE)
        {
            std::cerr<<"MATIO::load_to_arma(not a double matrix)"<<std::endl;
            return false;
        }
    }
    if(typeid(typename _mat::elem_type) == typeid(float))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_SINGLE,MAT_T_SINGLE,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(uint8_t))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT8,MAT_T_UINT8,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(long long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_INT64,MAT_T_INT64,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(unsigned long long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT64,MAT_T_UINT64,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(unsigned long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_UINT32,MAT_T_UINT32,2,dims,data,0);
    }
    if(typeid(typename _mat::elem_type) == typeid(long))
    {
        matvar = Mat_VarCreate(varname.c_str(),MAT_C_INT32,MAT_T_INT32,2,dims,data,0);
    }
    if(!matvar){
        Mat_Close(mat);
        std::cerr<<"MATIO::save_to_matlab(unknown type)"<<std::endl;
        return false;
    }
    Mat_VarWrite( mat, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(mat);
}
}
#endif // IOCORE_HPP
