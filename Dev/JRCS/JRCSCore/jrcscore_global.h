#ifndef JRCSCORE_GLOBAL_H
#define JRCSCORE_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(JRCSCORE_LIBRARY)
#  define JRCSCORESHARED_EXPORT Q_DECL_EXPORT
#else
#  define JRCSCORESHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // JRCSCORE_GLOBAL_H
