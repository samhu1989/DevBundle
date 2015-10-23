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

