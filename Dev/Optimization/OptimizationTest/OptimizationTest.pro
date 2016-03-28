#-------------------------------------------------
#
# Project created by QtCreator 2016-03-28T17:58:00
#
#-------------------------------------------------

QT       += testlib

QT       -= gui

TARGET = tst_optimizationtest
CONFIG   += console
CONFIG   -= app_bundle
DEFINES += USE_SSE
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread

TEMPLATE = app

SOURCES += tst_optimizationtest.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

DESTDIR = $$OUT_PWD/../../Dev_RunTime/bin

LIBS += -lopenblas

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

win32:CONFIG(release, debug|release): LIBS += -L$$DESTDIR/ -lOptimizationCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$DESTDIR/ -lOptimizationCore

INCLUDEPATH += $$PWD/../OptimizationCore
DEPENDPATH += $$PWD/../OptimizationCore
