TEMPLATE = app
TARGET = MolShaCS
QT += widgets
DEFINES += HAS_GUI
CONFIG += static
HEADERS += src/Mol2.h \
    src/GUI/QtWriter.h \
    src/CORREL.h \
    src/GUI/GUI.h \
    src/GUI/Widget.h \
    src/GUI/plotter.h \
    src/Grid.h \
    src/Minimizer2.h \
    src/Mol.h \
    src/Parser.h \
    src/RunEngine.h \
    src/Writer.h \
    src/Gaussian.h \
    src/main.h \
    src/ElSA.h
SOURCES += src/Mol2.cpp \
    src/GUI/QtWriter.cpp \
    src/CORREL.cpp \
    src/GUI/GUI.cpp \
    src/GUI/Widget.cpp \
    src/GUI/plotter.cpp \
    src/Grid.cpp \
    src/Minimizer2.cpp \
    src/Mol.cpp \
    src/Parser.cpp \
    src/RunEngine.cpp \
    src/Writer.cpp \
    src/Gaussian.cpp \
    src/main.cpp
RESOURCES += src/GUI/GUI.qrc
LIBS += -L/Users/asn/workspace/libela/trunk/lib/gsl/lib \
    -L/Users/asn/workspace/libela/trunk/lib/nlopt/lib \
    -L/Users/asn/workspace/libela/trunk/lib/zlib/lib \
    -lgsl -lgslcblas -lz -lnlopt_cxx
INCLUDEPATH += /Users/asn/workspace/libela/trunk/lib/gsl/include \
    /Users/asn/workspace/libela/trunk/lib/nlopt/include/ \
    /Users/asn/workspace/libela/trunk/lib/zlib/include/
QMAKE_CXXFLAGS = -O3 -ffast-math -std=c++11
CXXFLAGS=-std=c++11
