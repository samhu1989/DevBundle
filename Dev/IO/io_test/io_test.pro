#-------------------------------------------------
#
# Project created by QtCreator 2016-10-09T10:49:15
#
#-------------------------------------------------

QT       += testlib
QT       -= gui
TARGET = tst_io_testtest
CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread
TEMPLATE = app
SOURCES += \
    tst_io_test.cpp
DEFINES += SRCDIR=\\\"$$PWD/\\\"

DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

LIBS += -lopenblas

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib


win32: LIBS += -L$$DESTDIR/ -lIOCore

INCLUDEPATH += $$PWD/../IOCore
DEPENDPATH += $$PWD/../IOCore

LIBS += -lz -lhdf5
