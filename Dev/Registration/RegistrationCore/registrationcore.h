#ifndef REGISTRATIONCORE_H
#define REGISTRATIONCORE_H
#include "registrationcore_global.h"
#include <armadillo>
#include <iostream>
#include <memory>
#include "MeshType.h"
#include "RegistrationBase.h"
#include "jrmpc.h"
#include "coherentpointdrift.h"
#include "RegistrationThreadT.h"
//CPDRigid3D DefaultMesh Registration Thread:
template class REGISTRATIONCORESHARED_EXPORT Registration::CPDRigid3D<DefaultMesh>;
template class REGISTRATIONCORESHARED_EXPORT Registration::JRMPC<DefaultMesh>;
class REGISTRATIONCORESHARED_EXPORT CPDR3D_DM_R_Thread:public Registration::RegistrationThreadT<Registration::CPDRigid3D<DefaultMesh>,DefaultMesh>
{
    Q_OBJECT
public slots:
    void quit()
    {
        reg_.quit();
        Registration::RegistrationThreadT<Registration::CPDRigid3D<DefaultMesh>,DefaultMesh>::quit();
    }

protected:
    void run(void){compute();}
private:
};

class REGISTRATIONCORESHARED_EXPORT JRMPC_Thread:public Registration::RegistrationThreadT<Registration::JRMPC<DefaultMesh>,DefaultMesh>
{
    Q_OBJECT
public slots:
    void quit()
    {
        reg_.quit();
        Registration::RegistrationThreadT<Registration::JRMPC<DefaultMesh>,DefaultMesh>::quit();
    }
protected:
    void run(void){compute();}
private:
};
#endif // REGISTRATIONCORE_H
