#-------------------------------------------------
#
# Project created by QtCreator 2015-10-27T10:26:41
#
#-------------------------------------------------

QT       -= gui

TARGET = Common
TEMPLATE = lib

DEFINES += COMMON_LIBRARY
CONFIG += c++11
SOURCES += common.cpp \
    MeshColor.cpp

HEADERS += common.h\
        common_global.h \
    MeshColor.h \
    MeshColor.hpp \
    MeshType.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DESTDIR = $$OUT_PWD/../../Dev_RunTime/bin

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll
