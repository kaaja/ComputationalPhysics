TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    jacobi.cpp \
    project2.cpp \  #\
    eigenvalueBisection.cpp \
    lanczos.cpp
#    tests-main.cpp \
#    test-functions.cpp

INCLUDEPATH += /usr/local/include
LIBS += -L/usr/local/lib
LIBS += -larmadillo -llapack -lblas

HEADERS += \
    jacobi.h \
    catch.hpp \
    eigenvalueBisection.h \
    lanczos.h
