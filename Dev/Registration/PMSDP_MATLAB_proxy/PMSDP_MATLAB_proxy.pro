#-------------------------------------------------
#
# Project created by QtCreator 2016-07-19T19:31:26
#
#-------------------------------------------------

QT       -= gui

TARGET = PMSDP_MATLAB_proxy
TEMPLATE = lib

DEFINES += PMSDP_MATLAB_PROXY_LIBRARY

DEFINES += __MW_STDINT_H__

SOURCES += pmsdp_matlab_proxy.cpp

HEADERS += pmsdp_matlab_proxy.h\
        pmsdp_matlab_proxy_global.h

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

win32: LIBS += -L$$PWD/../../../3rdParty/MATLAB/extern/lib/win64/microsoft -llibmx -lmclmcr -lmclmcrrt

INCLUDEPATH += $$PWD/../../../3rdParty/MATLAB/extern/include
DEPENDPATH += $$PWD/../../../3rdParty/MATLAB/extern/include

win32: LIBS += -L$$PWD/../../../3rdParty/MATLAB/lib/ -llibPMSDP

INCLUDEPATH += $$PWD/../../../3rdParty/MATLAB/include
DEPENDPATH += $$PWD/../../../3rdParty/MATLAB/include
