#ifndef FEATURECORE_H
#define FEATURECORE_H
#include <memory>
#include "featurecore_global.h"
#include "pointnormal.h"
namespace  Feature
{
template void FEATURECORESHARED_EXPORT computeNormalPoint<DefaultMesh>(
                                                                DefaultMesh& mesh,
                                                                float r,
                                                                int k
                                                                        );
}

#endif // FEATURECORE_H
