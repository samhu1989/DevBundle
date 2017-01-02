#-------------------------------------------------
#
# Project created by QtCreator 2017-01-02T16:44:06
#
#-------------------------------------------------

QT       -= gui

TARGET = DL
TEMPLATE = lib

DEFINES += DL_LIBRARY

SOURCES += \
    net.cpp \
    layer.cpp \
    linearlayer.cpp \
    relulayer.cpp

HEADERS += dl.h\
        dl_global.h \
    net.h \
    layer.h \
    linearlayer.h \
    relulayer.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

win32: LIBS += -L$$DESTDIR/ -lCommon

INCLUDEPATH += $$PWD/../../Common
DEPENDPATH += $$PWD/../../Common

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../3rdParty/SuperLU/include

win32: LIBS += -L$$DESTDIR/ -lOptimizationCore

INCLUDEPATH += $$PWD/../../Optimization/OptimizationCore
DEPENDPATH += $$PWD/../../Optimization/OptimizationCore

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

LIBS += -lopenblas
