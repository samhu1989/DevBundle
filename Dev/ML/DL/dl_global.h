#ifndef DL_GLOBAL_H
#define DL_GLOBAL_H
#include <QtCore/qglobal.h>
#if defined(DL_LIBRARY)
#  define DLSHARED_EXPORT Q_DECL_EXPORT
#else
#  define DLSHARED_EXPORT Q_DECL_IMPORT
#endif

#endif // DL_GLOBAL_H
