TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    isingmodel.cpp \
    lib.cpp

HEADERS += \
    isingmodel.h \
    lib.h
