#ifndef IOCORE_HPP
#define IOCORE_HPP
#include "iocore.h"
#include "matio.h"
#include <typeinfo>
#include <memory>
#include <assert.h>
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
void save_to_matlab(const std::vector<std::shared_ptr<_mat>>& m_lst, const std::string &file, const std::vector<std::string> &var)
{
    assert(m_lst.size()==var.size());
    mat_t *mat;
    mat = Mat_CreateVer(file.c_str(),NULL,MAT_FT_MAT73);
    for(int i = 0 ; i < m_lst.size() ; ++i )
    {
        const arma::mat& m = (*m_lst[i]);
        void* data = (void*)m.memptr();
        size_t dims[2] = {m.n_rows,m.n_cols};
        matvar_t *matvar = NULL;
        std::string varname = var[i];
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
    }
    Mat_Close(mat);
}

template<typename _mat>
bool load_to_arma(_mat& m,const std::string& file,const std::string& var)
{
//    std::cerr<<"load_to_arma:a"<<std::endl;
    mat_t *mat = Mat_Open(file.c_str(),MAT_ACC_RDONLY);
    if(!mat){
        std::cerr<<"MATIO::load_to_arma(failed to open "<<file<<")"<<std::endl;
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
    matvar = Mat_VarRead(mat,varname.c_str());
//    std::cerr<<"load_to_arma:b"<<std::endl;
    if(!matvar)
    {
        std::cerr<<"MATIO::load_to_arma(failed to find variable "<<varname<<")"<<std::endl;
        Mat_Close(mat);
        return false;
    }
    std::shared_ptr<char> t_name;
    const char* s_name;
    switch(matvar->data_type)
    {
    case MAT_T_DOUBLE:
        s_name = typeid(double).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_SINGLE:
        s_name = typeid(float).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_UINT64:
        s_name = typeid(uint64_t).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_INT64:
        s_name = typeid(int64_t).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_UINT32:
        s_name = typeid(uint32_t).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_INT32:
        s_name = typeid(int32_t).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    case MAT_T_UINT8:
        s_name = typeid(uint8_t).name();
        t_name.reset(new char[std::strlen(s_name)+1]);
        std::strcpy(t_name.get(),s_name);
        break;
    default:
        std::cerr<<"MATIO::load_to_arma(unknown data_type)"<<std::endl;
        Mat_VarFree(matvar);
        Mat_Close(mat);
        return false;
    }
//    std::cerr<<"load_to_arma:c"<<std::endl;
//    std::cerr<<"l0:"<<std::strlen(t_name.get())<<std::endl;
    if(0==std::strcmp(typeid(typename _mat::elem_type).name(),t_name.get()))
    {
//        std::cerr<<"load_to_arma:c2"<<std::endl;
        if( matvar->rank != 2 )
        {
            std::cerr<<"MATIO::load_to_arma(not a matrix,a data with rank="<<matvar->rank<<")"<<std::endl;
            Mat_VarFree(matvar);
            Mat_Close(mat);
            return false;
        }
        if( matvar->dims[0] != 1 && m.is_row)
        {
            std::cerr<<"MATIO::load_to_arma(expect a row vec but get a dims[0]="<<matvar->dims[0]<<")"<<std::endl;
            Mat_VarFree(matvar);
            Mat_Close(mat);
            return false;
        }else if( matvar->dims[1] != 1 && m.is_col )
        {
            std::cerr<<"MATIO::load_to_arma(expect a col vec but get a dims[1]="<<matvar->dims[1]<<")"<<std::endl;
            Mat_VarFree(matvar);
            Mat_Close(mat);
            return false;
        }else{
            if(m.is_row)m = _mat(1,matvar->dims[1],arma::fill::zeros);
            else if(m.is_col)m = _mat(matvar->dims[0],1,arma::fill::zeros);
            else m = _mat(matvar->dims[0],matvar->dims[1],arma::fill::zeros);
            std::memcpy((void*)m.memptr(),matvar->data,matvar->data_size*matvar->dims[0]*matvar->dims[1]);
        }
    }else{
        std::cerr<<"MATIO::load_to_arma(not a "<<typeid(typename _mat::elem_type).name()<<" matrix)"<<std::endl;
        Mat_VarFree(matvar);
        Mat_Close(mat);
        return false;
    }
    Mat_VarFree(matvar);
    Mat_Close(mat);
    return true;
}

}
#endif // IOCORE_HPP
