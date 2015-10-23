#-------------------------------------------------
#
# Project created by QtCreator 2015-10-21T16:28:22
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = OpenMeshViewer
TEMPLATE = app
CONFIG += c++11
DEFINES += OM_STATIC_BUILD
DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

SOURCES +=\
    meshviewer.cc

HEADERS  +=

FORMS    +=

win32:CONFIG(release, debug|release): LIBS += -lfreeglut
else:win32:CONFIG(debug, debug|release): LIBS += -lfreeglut

win32:CONFIG(release, debug|release): LIBS += -L$$OUT_PWD/../VisualizationCore/release/ -lVisualizationCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$OUT_PWD/../VisualizationCore/debug/ -lVisualizationCore

INCLUDEPATH += $$PWD/../VisualizationCore
DEPENDPATH += $$PWD/../VisualizationCore

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCore.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCored.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCore.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCored.lib

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshTools
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshTools.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshToolsd.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshTools.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshToolsd.lib
