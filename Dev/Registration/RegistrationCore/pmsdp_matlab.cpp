#include "pmsdp_matlab.h"
#include <QDebug>
#include <QFileInfo>
#include <QCoreApplication>
namespace Registration{
PMSDP_MATLAB::PMSDP_MATLAB(void):_compute(NULL),_init(NULL),_terminate(NULL)
{
    QCoreApplication::addLibraryPath("./bin/");
    QCoreApplication::addLibraryPath("D:/MATLAB/R2014/bin/win64/");
    QCoreApplication::addLibraryPath("D:/MATLAB/R2014/runtime/win64/");
    lib.setFileName("./bin/PMSDP_MATLAB_proxy.dll");
    if(!lib.load())
    {
        qDebug()<<"Failed to Load "<<lib.fileName();
        qDebug()<<"With Error:"<<lib.errorString();
        lib.setFileName("./PMSDP_MATLAB_proxy.dll");
    }
    if(lib.load()){
        _init = (InitFunc)lib.resolve("initMatlab");
        _terminate = (EndFunc)lib.resolve("terminateMatlab");
        _compute = (Func)lib.resolve("compute");
    }
    else{
        qDebug()<<"Failed to Load "<<lib.fileName();
        qDebug()<<"With Error:"<<lib.errorString();
    }
}
PMSDP_MATLAB::PMSDP_MATLAB(const QString& fileName):_compute(NULL),_init(NULL),_terminate(NULL),lib(fileName)
{
    QCoreApplication::addLibraryPath("./bin/");
    QCoreApplication::addLibraryPath("D:/MATLAB/R2014/bin/win64/");
    QCoreApplication::addLibraryPath("D:/MATLAB/R2014/runtime/win64/");
    QFileInfo info(fileName);
    if(!info.exists())
    {
        qDebug()<<info.absoluteFilePath()<<" doesn't exist";
    }
    if(!lib.load())
    {
        qDebug()<<"Failed to Load "<<lib.fileName();
        qDebug()<<"With Error:"<<lib.errorString();
        lib.setFileName("./bin/PMSDP_MATLAB_proxy.dll");
        if(!lib.load())
        {
            qDebug()<<"Failed to Load "<<lib.fileName();
            qDebug()<<"With Error:"<<lib.errorString();
            lib.setFileName("./PMSDP_MATLAB_proxy.dll");
        }
    }
    if(lib.isLoaded()){
        _init = (InitFunc)lib.resolve("initMatlab");
        _terminate = (EndFunc)lib.resolve("terminateMatlab");
        _compute = (Func)lib.resolve("compute");
    }
    else {
        qDebug()<<lib.fileName()<<" is still not Loaded";
    }
}

bool PMSDP_MATLAB::init()
{
    if(!lib.isLoaded())return false;
    if(!_init)return false;
    return (*_init)();
}

void PMSDP_MATLAB::compute(
        const arma::mat& P,
        const arma::mat& Q,
        arma::mat& R,
        arma::uvec& X
        )
{
    if(NULL==_compute)
    {
        std::cerr<<"Failed to resolve \"compute\" from "<<lib.fileName().toStdString()<<std::endl;
        qDebug()<<"With Error:"<<lib.errorString();
        return;
    }
    (*_compute)(P,Q,R,X);
}
}
