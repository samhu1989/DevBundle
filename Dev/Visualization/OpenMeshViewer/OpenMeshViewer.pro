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


win32:CONFIG(release, debug|release): LIBS += -L$$DESTDIR/ -lVisualizationCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$DESTDIR/ -lVisualizationCore

INCLUDEPATH += $$PWD/../VisualizationCore
DEPENDPATH += $$PWD/../VisualizationCore

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll

win32:CONFIG(release, debug|release): LIBS += -lfreeglut
else:win32:CONFIG(debug, debug|release): LIBS += -lfreeglut
