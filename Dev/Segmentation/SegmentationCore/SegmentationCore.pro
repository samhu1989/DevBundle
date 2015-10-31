#-------------------------------------------------
#
# Project created by QtCreator 2015-10-31T19:07:30
#
#-------------------------------------------------

QT       -= gui

TARGET = SegmentationCore
TEMPLATE = lib

DEFINES += SEGMENTATIONCORE_LIBRARY

SOURCES += segmentationcore.cpp

HEADERS += segmentationcore.h\
        segmentationcore_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
