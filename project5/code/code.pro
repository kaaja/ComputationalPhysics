TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp \
    twodimensionaldiffusionsolver.cpp \
    lib.cpp

HEADERS += \
    solver.h \
    twodimensionaldiffusionsolver.h \
    lib.h

LIBS += -larmadillo -llapack -lblas

QMAKE_CXXFLAGS+= -fopenmp
QMAKE_LFLAGS +=  -fopenmp
