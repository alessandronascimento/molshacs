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
LIBS += C:/Users/Nascimento/workspace/MolShaCS_2/lib/gsl/lib/libgsl.a \
    C:/Users/Nascimento/workspace/MolShaCS_2/lib/gsl/lib/libgslcblas.a \
    C:/Users/Nascimento/workspace/MolShaCS_2/lib/nlopt/lib/libnlopt.a \
    C:/Users/Nascimento/workspace/MolShaCS_2/lib/zlib/lib/libz.a
INCLUDEPATH += C:/Users/Nascimento/workspace/MolShaCS_2/lib/gsl/include \
    C:/Users/Nascimento/workspace/MolShaCS_2/lib/nlopt/include/ \
    C:/Users/Nascimento/workspace/MolShaCS_2/lib/zlib/include/
QMAKE_CXXFLAGS = -O3 -ffast-math -std=c++11 -static
CXXFLAGS=-std=c++11
