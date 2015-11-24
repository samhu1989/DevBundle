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
#DEFINES += OM_STATIC_BUILD
DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

SOURCES += \
    meshviewer.cpp \
    mainwindow.cpp \
    computenormalthread.cpp \
    computeoctreethread.cpp \
    computesupervoxelthread.cpp

HEADERS  += \
    mainwindow.h \
    computenormalthread.h \
    computeoctreethread.h \
    computesupervoxelthread.h

FORMS    += \
    mainwindow.ui


win32: LIBS += -L$$DESTDIR/ -lFeatureCore

INCLUDEPATH += $$PWD/../../Feature/FeatureCore
DEPENDPATH += $$PWD/../../Feature/FeatureCore

win32:CONFIG(release, debug|release): LIBS += -L$$DESTDIR/ -lVisualizationCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$DESTDIR/ -lVisualizationCore

INCLUDEPATH += $$PWD/../VisualizationCore
DEPENDPATH += $$PWD/../VisualizationCore

win32: LIBS += -L$$DESTDIR/ -lCommon

INCLUDEPATH += $$PWD/../../Common
DEPENDPATH += $$PWD/../../Common

win32:CONFIG(release, debug|release): LIBS += -lfreeglut
else:win32:CONFIG(debug, debug|release): LIBS += -lfreeglut

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

INCLUDEPATH += $$PWD/../../../3rdParty/NanoFlann/include
DEPENDPATH += $$PWD/../../../3rdParty/NanoFlann/include

LIBS += -lopenblas

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lblas
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lblas

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libblas.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libblas.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/blas.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/blas.lib



win32: LIBS += -L$$DESTDIR/ -lSegmentationCore

INCLUDEPATH += $$PWD/../../Segmentation/SegmentationCore
DEPENDPATH += $$PWD/../../Segmentation/SegmentationCore
