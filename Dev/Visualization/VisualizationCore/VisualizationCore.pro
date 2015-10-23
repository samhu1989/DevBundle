#-------------------------------------------------
#
# Project created by QtCreator 2015-10-23T16:39:51
#
#-------------------------------------------------

QT       += widgets opengl

TARGET = VisualizationCore
TEMPLATE = lib

DEFINES += VISUALIZATIONCORE_LIBRARY
DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin
SOURCES += visualizationcore.cpp \
    MeshViewerWidgetT.cc \
    QGLViewerWidget.cc

HEADERS += visualizationcore.h\
        visualizationcore_global.h \
    MeshViewerWidget.hh \
    MeshViewerWidgetT.hh \
    QGLViewerWidget.hh

unix {
    target.path = /usr/lib
    INSTALLS += target
}

win32:CONFIG(release, debug|release): LIBS += -lfreeglut
else:win32:CONFIG(debug, debug|release): LIBS += -lfreeglut

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

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
