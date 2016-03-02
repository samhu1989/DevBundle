#ifndef REGISTRATIONBASE
#define REGISTRATIONBASE
#include "MeshType.h"
#include "registrationcore_global.h"
namespace Registration
{
    class REGISTRATIONCORESHARED_EXPORT RegistrationBase
    {
    public:
        RegistrationBase():end_(false){}
        virtual bool initForThread(void*){return false;}
        virtual void compute()=0;
        virtual const std::string& errorForThread(void)const {return error_string_;}
        virtual void quit(){end_=true;}
        static void alignAroundZ(const arma::fmat&,const arma::fmat&,arma::fmat&);
        static void alignAroundZ(
                const arma::fmat&,
                const arma::fmat&,
                const arma::fmat&,
                const arma::fmat&,
                arma::fmat&
                );
        static double closestError(const arma::fmat&,const arma::fmat&);
        static double closestError(
                const arma::fmat&,
                const arma::fmat&,
                const arma::fmat&,
                const arma::fmat&
                );
    protected:
        std::string error_string_;
        bool end_;
    };
}
#endif // REGISTRATIONBASE

