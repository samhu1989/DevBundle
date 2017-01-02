#-------------------------------------------------
#
# Project created by QtCreator 2016-03-28T16:20:47
#
#-------------------------------------------------

QT       -= gui

TARGET = OptimizationCore
TEMPLATE = lib

DEFINES += OPTIMIZATIONCORE_LIBRARY
DEFINES += USE_SSE HAVE_CONFIG_H
DEFINES -= __SSE3__
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread
SOURCES += optimizationcore.cpp \
    lbfgs.cpp \
    LBFGS/lbfgscore.c \
    sdp.cpp

HEADERS += optimizationcore.h\
        optimizationcore_global.h \
    lbfgs.h \
    LBFGS/arithmetic_ansi.h \
    LBFGS/arithmetic_sse_double.h \
    LBFGS/arithmetic_sse_float.h \
    LBFGS/lbfgscore.h \
    sdp.h \
    LBFGS/config.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib

INCLUDEPATH += $$PWD/../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../3rdParty/SuperLU/include

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libblas.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libblas.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/blas.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/blas.lib

win32: LIBS += -L$$PWD/../../../3rdParty/CSDP/lib/ -lsdp

INCLUDEPATH += $$PWD/../../../3rdParty/CSDP/include
DEPENDPATH += $$PWD/../../../3rdParty/CSDP/include

win32:!win32-g++: PRE_TARGETDEPS += $$PWD/../../../3rdParty/CSDP/lib/sdp.lib
else:win32-g++: PRE_TARGETDEPS += $$PWD/../../../3rdParty/CSDP/lib/libsdp.a

LIBS += -lopenblas -lm
