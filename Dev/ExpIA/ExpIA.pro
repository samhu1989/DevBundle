#-------------------------------------------------
#
# Project created by QtCreator 2015-11-11T21:32:47
#
#-------------------------------------------------

QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ExpIA
TEMPLATE = app
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread

SOURCES += main.cpp\
        mainwindow.cpp \
    regiongrowthread.cpp \
    unifylabelcolorsizethread.cpp \
    unifylabelmannual.cpp \
    patchpairview.cpp \
    updateobjectmodel.cpp \
    objectmodel.cpp

HEADERS  += mainwindow.h \
    regiongrowthread.h \
    unifylabelcolorsizethread.h \
    unifylabelmannual.h \
    patchpairview.h \
    updateobjectmodel.h \
    objectmodel.h

FORMS    += mainwindow.ui \
    unifylabelmannual.ui \
    patchpairview.ui \
    updateobjectmodel.ui

DESTDIR = $$OUT_PWD/../../Dev_RunTime/bin

win32: LIBS += -L$$DESTDIR/ -lCommon

INCLUDEPATH += $$PWD/../Common
DEPENDPATH += $$PWD/../Common

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshCore.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshCored.dll

INCLUDEPATH += $$PWD/../../3rdParty/OpenMesh/include
DEPENDPATH += $$PWD/../../3rdParty/OpenMesh/include

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshTools.dll
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/OpenMesh/lib/ -lOpenMeshToolsd.dll

LIBS += -lopenblas

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/superlu.lib


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../3rdParty/SuperLU/lib/ -lblas
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../3rdParty/SuperLU/lib/ -lblas

INCLUDEPATH += $$PWD/../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libblas.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/libblas.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/blas.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../3rdParty/SuperLU/lib/blas.lib

win32: LIBS += -L$$DESTDIR/ -lVisualizationCore

INCLUDEPATH += $$PWD/../Visualization/VisualizationCore
DEPENDPATH += $$PWD/../Visualization/VisualizationCore

win32: LIBS += -L$$DESTDIR/ -lSegmentationCore

INCLUDEPATH += $$PWD/../Segmentation/SegmentationCore
DEPENDPATH += $$PWD/../Segmentation/SegmentationCore

INCLUDEPATH += $$PWD/../../3rdParty/NanoFlann/include
DEPENDPATH += $$PWD/../../3rdParty/NanoFlann/include

win32: LIBS += -L$$DESTDIR/ -lFeatureCore

INCLUDEPATH += $$PWD/../Feature/FeatureCore
DEPENDPATH += $$PWD/../Feature/FeatureCore

win32: LIBS += -L$$DESTDIR/ -lRegistrationCore

INCLUDEPATH += $$PWD/../Registration/RegistrationCore
DEPENDPATH += $$PWD/../Registration/RegistrationCore
