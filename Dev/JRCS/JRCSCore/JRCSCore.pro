#-------------------------------------------------
#
# Project created by QtCreator 2016-04-03T15:36:15
#
#-------------------------------------------------

QT  += core  gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = JRCSCore
TEMPLATE = lib
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread
DEFINES += JRCSCORE_LIBRARY
#DEFINES -= UNICODE
#DEFINES -= _UNICODE
SOURCES += \
    jrcsbase.cpp \
    jrcsinitbase.cpp \
    jrcsinitexternal.cpp \
    jrcsaoni.cpp \
    jrcsaopt.cpp \
    sjrcsbase.cpp \
    jrcsbilateral.cpp \
    jrcsprimitive.cpp \
    jrcscube.cpp \
    jrcsbox.cpp

HEADERS +=\
        jrcscore_global.h \
    jrcsbase.h \
    jrcsinitbase.h \
    jrcsinitexternal.h \
    jrcsaoni.h \
    jrcsaopt.h \
    sjrcsbase.h \
    jrcsbilateral.h \
    jrcsprimitive.h \
    jrcscube.h \
    jrcsbox.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

win32: LIBS += -L$$DESTDIR/ -lCommon

INCLUDEPATH += $$PWD/../../Common
DEPENDPATH += $$PWD/../../Common

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCore
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/OpenMesh/lib/ -lOpenMeshCored

INCLUDEPATH += $$PWD/../../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../../3rdParty/OpenMesh/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCore.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCored.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCore.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCored.lib

LIBS += -lopenblas

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib


#win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lblas
#else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lblas

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

#win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libblas.a
#else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libblas.a
#else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/blas.lib
#else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/blas.lib

INCLUDEPATH += $$PWD/../../../3rdParty/NanoFlann/include
DEPENDPATH += $$PWD/../../../3rdParty/NanoFlann/include

win32: LIBS += -L$$OUT_PWD/../../Segmentation/SegmentationCore/ -lSegmentationCore

INCLUDEPATH += $$PWD/../../Segmentation/SegmentationCore
DEPENDPATH += $$PWD/../../Segmentation/SegmentationCore

win32: LIBS += -L$$OUT_PWD/../../Feature/FeatureCore/ -lFeatureCore

INCLUDEPATH += $$PWD/../../Feature/FeatureCore
DEPENDPATH += $$PWD/../../Feature/FeatureCore

win32: LIBS += -L$$OUT_PWD/../../ML/Clustering/ -lClustering

INCLUDEPATH += $$PWD/../../ML/Clustering
DEPENDPATH += $$PWD/../../ML/Clustering

win32: LIBS += -L$$OUT_PWD/../../IO/IOCore/ -lIOCore

INCLUDEPATH += $$PWD/../../IO/IOCore
DEPENDPATH += $$PWD/../../IO/IOCore
