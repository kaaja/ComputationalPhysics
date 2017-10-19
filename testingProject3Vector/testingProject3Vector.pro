TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    ../project3/planetmercury.cpp \
    ../project3/planet.cpp \
    ../project3/solver.cpp

HEADERS += \
    planet.h \
    planetmercury.h \
    solver.h
