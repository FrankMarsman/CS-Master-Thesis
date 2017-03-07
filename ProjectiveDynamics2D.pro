#-------------------------------------------------
#
# Project created by QtCreator 2017-03-07T11:17:48
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

CONFIG += c++11

TARGET = ProjectiveDynamics2D
TEMPLATE = app


QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
QMAKE_CFLAGS_RELEASE += -fopenmp

INCLUDEPATH += "D:/Thesis TUD/Libraries/Eigen/eigen-eigen-26667be4f70b/eigen-eigen-26667be4f70b/Eigen"

SOURCES += main.cpp\
        mainwindow.cpp \
    auxvar2d.cpp \
    sim2d.cpp \
    qcustomplot.cpp

HEADERS  += mainwindow.h \
    auxvar2d.h \
    sim2d.h \
    qcustomplot.h

FORMS    += mainwindow.ui
