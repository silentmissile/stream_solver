#-------------------------------------------------
#
# Project created by QtCreator 2013-03-06T15:51:17
#
#-------------------------------------------------

INCLUDEPATH += ~/Develop/eigen-eigen-3.1.2/
QT       -= gui

TARGET = stream_solver
TEMPLATE = lib

DEFINES += STREAM_SOLVER_LIBRARY

SOURCES += stream_solver.cpp \
    spline.cpp

HEADERS += stream_solver.h\
        stream_solver_global.h \
    spline.h

unix:!symbian {
    maemo5 {
        target.path = /opt/usr/lib
    } else {
        target.path = /usr/lib
    }
    INSTALLS += target
}
