#-------------------------------------------------
#
# Project created by QtCreator 2018-04-27T17:15:58
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = ManipQT
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SUBDIRS += tests
SOURCES += \
        main.cpp \
        mainwindow.cpp \
    alglib/alglibinternal.cpp \
    alglib/alglibmisc.cpp \
    alglib/ap.cpp \
    alglib/dataanalysis.cpp \
    alglib/diffequations.cpp \
    alglib/fasttransforms.cpp \
    alglib/integration.cpp \
    alglib/interpolation.cpp \
    alglib/linalg.cpp \
    alglib/optimization.cpp \
    alglib/solvers.cpp \
    alglib/specialfunctions.cpp \
    alglib/statistics.cpp \
    manipcalculations.cpp \
    data_saver.cpp \
    omp_settings.cpp \
    about.cpp \
    splain_calculations.cpp

HEADERS += \
        mainwindow.h \
    alglib/alglibinternal.h \
    alglib/alglibmisc.h \
    alglib/ap.h \
    alglib/dataanalysis.h \
    alglib/diffequations.h \
    alglib/fasttransforms.h \
    alglib/integration.h \
    alglib/interpolation.h \
    alglib/linalg.h \
    alglib/optimization.h \
    alglib/solvers.h \
    alglib/specialfunctions.h \
    alglib/statistics.h \
    alglib/stdafx.h \
    manipcalculations.h \
    data_saver.h \
    omp_settings.h \
    about.h \
    splain_calculations.h

FORMS += \
        mainwindow.ui \
    omp_settings.ui \
    about.ui

RESOURCES += \
    resourses.qrc
#need add in QMake build arguments "QMAKE_LIBS+=-static -lgomp -lpthread" "QMAKE_CXXFLAGS+=-fopenmp"
QMAKE_LIBS+=-static -lgomp -lpthread
QMAKE_CXXFLAGS+=-fopenmp
QMAKE_LFLAGS += -fopenmp
QMAKE_CFLAGS_DEBUG += -fopenmp
#QMAKE_CFLAGS_RELEASE += -fopenmp
