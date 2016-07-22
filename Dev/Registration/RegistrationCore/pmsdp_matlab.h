#ifndef PMSDP_MATLAB_H
#define PMSDP_MATLAB_H
#include <QLibrary>
#include <armadillo>
#include <QDebug>
#include <assert.h>
namespace Registration{
class PMSDP_MATLAB
{
public:
    typedef void (*Func)(
            const arma::mat &P,
            const arma::mat &Q,
            arma::mat& R,
            arma::uvec& X
            );
    typedef bool (*InitFunc)(void);
    typedef void (*EndFunc)(void);
    PMSDP_MATLAB(void);
    PMSDP_MATLAB(const QString& fileName);
    bool init(void);
    void compute(const arma::mat &P,
            const arma::mat &Q,
            arma::mat& R,
            arma::uvec& X
            );
    ~PMSDP_MATLAB(){
        if(lib.isLoaded())
        {
            assert(_terminate);
        }
        if(_terminate);
        {
            (*_terminate)();
        }
        if(!lib.unload())
        {
            qDebug()<<lib.errorString();
        }
    }
private:
    QLibrary lib;
    InitFunc _init;
    Func _compute;
    EndFunc _terminate;
};
}
#endif // PMSDP_MATLAB_H
