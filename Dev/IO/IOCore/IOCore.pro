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
        iocore_global.h \
    iocore.hpp

unix {
    target.path = /usr/lib
    INSTALLS += target
}

DESTDIR = $$OUT_PWD/../../../Dev_RunTime/bin

LIBS += -lopenblas

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../3rdParty/SuperLU/lib/ -lsuperlu

INCLUDEPATH += $$PWD/../../../3rdParty/SuperLU/include
DEPENDPATH += $$PWD/../../../3rdParty/SuperLU/include

win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/libsuperlu.a
else:win32:!win32-g++:CONFIG(release, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib
else:win32:!win32-g++:CONFIG(debug, debug|release): PRE_TARGETDEPS += $$PWD/../../../3rdParty/SuperLU/lib/superlu.lib

win32: LIBS += -L$$PWD/../../../3rdParty/MATIO/lib/ -llibmatio

INCLUDEPATH += $$PWD/../../../3rdParty/MATIO/include
DEPENDPATH += $$PWD/../../../3rdParty/MATIO/include

LIBS += -lz -lhdf5
