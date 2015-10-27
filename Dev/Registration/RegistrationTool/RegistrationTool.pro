QT  += core gui opengl
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = RegistrationTool
CONFIG += console
CONFIG += c++11
DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin


TEMPLATE = app

SOURCES += main.cpp \
    mainwindow.cpp

FORMS += \
    mainwindow.ui

HEADERS += \
    mainwindow.h


win32: LIBS += -L$$DESTDIR/ -lVisualizationCore

INCLUDEPATH += $$PWD/../../Visualization/VisualizationCore
DEPENDPATH += $$PWD/../../Visualization/VisualizationCore

win32: LIBS += -L$$DESTDIR/ -lCommon

INCLUDEPATH += $$PWD/../../Common
DEPENDPATH += $$PWD/../../Common

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll

win32:CONFIG(release, debug|release): LIBS += -lfreeglut
else:win32:CONFIG(debug, debug|release): LIBS += -lfreeglut
