#-------------------------------------------------
#
# Project created by QtCreator 2016-10-09T14:59:42
#
#-------------------------------------------------

QT       -= gui

TARGET = IOCore
TEMPLATE = lib

DEFINES += IOCORE_LIBRARY

SOURCES += iocore.cpp

HEADERS += iocore.h\
        iocore_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
