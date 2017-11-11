TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    solver.cpp

HEADERS += \
    solver.h

LIBS += -larmadillo -llapack -lblas
