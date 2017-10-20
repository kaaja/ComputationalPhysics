TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    solver.cpp \
    planetmercury.cpp \
    alternativeplanet.cpp

HEADERS += \
    planet.h \
    solver.h \
    planetmercury.h \
    alternativeplanet.h
