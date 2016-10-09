#ifndef IOCORE_GLOBAL_H
#define IOCORE_GLOBAL_H

#include <QtCore/qglobal.h>

#if defined(IOCORE_LIBRARY)
#  define IOCORESHARED_EXPORT Q_DECL_EXPORT
#else
#  define IOCORESHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // IOCORE_GLOBAL_H
