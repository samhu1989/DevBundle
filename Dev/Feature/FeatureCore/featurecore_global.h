#ifndef FEATURECORE_GLOBAL_H
#define FEATURECORE_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(FEATURECORE_LIBRARY)
#  define FEATURECORESHARED_EXPORT Q_DECL_EXPORT
#else
#  define FEATURECORESHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // FEATURECORE_GLOBAL_H
