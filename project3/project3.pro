TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    solver.cpp \
    system.cpp \
    planetgeneralrelativityforce.cpp \
    planetalternativeforce.cpp

HEADERS += \
    planet.h \
    solver.h \
    system.h \
    planetgeneralrelativityforce.h \
    planetalternativeforce.h
