TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    twodimensionaldiffusionsolver.cpp

HEADERS += \
    solver.h \
    twodimensionaldiffusionsolver.h

LIBS += -larmadillo -llapack -lblas
