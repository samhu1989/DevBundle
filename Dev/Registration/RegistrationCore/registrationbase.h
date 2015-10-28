#ifndef REGISTRATIONBASE
#define REGISTRATIONBASE
#include "MeshType.h"
#include "registrationcore_global.h"
namespace Registration
{
    class REGISTRATIONCORESHARED_EXPORT RegistrationBase
    {
    public:
        virtual bool initForThread(void*)=0;
        virtual void compute()=0;
        virtual const std::string& errorForThread(void)const {return error_string_;}
    protected:
        std::string error_string_;
    };
}
#endif // REGISTRATIONBASE

