#-------------------------------------------------
#
# Project created by QtCreator 2015-10-31T19:07:30
#
#-------------------------------------------------

QT       += core gui

TARGET = SegmentationCore
TEMPLATE = lib
CONFIG += c++11
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp -lpthread
DEFINES += SEGMENTATIONCORE_LIBRARY
DEFINES += USE_64_BIT_PTR_CAST
SOURCES += segmentationcore.cpp \
    graphcut.cpp \
    MRF/src/BP-S.cpp \
    MRF/src/GCoptimization.cpp \
    MRF/src/graph.cpp \
    MRF/src/ICM.cpp \
    MRF/src/LinkedBlockList.cpp \
    MRF/src/maxflow.cpp \
    MRF/src/MaxProdBP.cpp \
    MRF/src/mrf.cpp \
    MRF/src/regions-maxprod.cpp \
    MRF/src/TRW-S.cpp \
    SAC/sac_plane.cpp \
    sac_segmentation.cpp \
    SAC/sac_parallel_plane.cpp \
    SAC/sac_perpendicular_plane.cpp \
    CRF/labelcompatibility.cpp \
    CRF/objective.cpp \
    CRF/pairwise.cpp \
    CRF/permutohedral.cpp \
    CRF/unary.cpp \
    CRF/util.cpp \
    densecrf.cpp \
    densecrf3d.cpp \
    hierarchicalization.cpp \
    pcaplaneequ.cpp

HEADERS += segmentationcore.h\
        segmentationcore_global.h \
    supervoxelclustering.h \
    supervoxelclustering.hpp \
    segmentationbase.h \
    regiongrowing.h \
    regiongrowing.hpp \
    graphcut.h \
    MRF/include/block.h \
    MRF/include/BP-S.h \
    MRF/include/energy.h \
    MRF/include/GCoptimization.h \
    MRF/include/graph.h \
    MRF/include/ICM.h \
    MRF/include/LinkedBlockList.h \
    MRF/include/MaxProdBP.h \
    MRF/include/mrf.h \
    MRF/include/regions-new.h \
    MRF/include/TRW-S.h \
    MRF/include/typeTruncatedQuadratic2D.h \
    sac_segmentation.h \
    SAC/sac_model.h \
    SAC/sac_plane.h \
    SAC/sac_parallel_plane.h \
    SAC/sac_perpendicular_plane.h \
    CRF/labelcompatibility.h \
    CRF/objective.h \
    CRF/pairwise.h \
    CRF/permutohedral.h \
    CRF/unary.h \
    CRF/util.h \
    densecrf.h \
    densecrf3d.h \
    hierarchicalization.h \
    pcaplaneequ.h \
    regiongrowingrgb.h \
    regiongrowingrgb.hpp \
    normalizedcuts.h \
    normalizedcuts.hpp \
    robustregion.hpp \
    robustregion.h

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

INCLUDEPATH += $$PWD/../../../3rdParty/NanoFlann/include
DEPENDPATH += $$PWD/../../../3rdParty/NanoFlann/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCore.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/libOpenMeshCored.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCore.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/OpenMesh/lib/OpenMeshCored.lib

LIBS += -lopenblas -larpack

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


win32: LIBS += -L$$DESTDIR/ -lClustering

INCLUDEPATH += $$PWD/../../ML/Clustering
DEPENDPATH += $$PWD/../../ML/Clustering
