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
RESOURCES += src/GUI/GUI.qrc \
    src/GUI/GUI.qrc
LIBS += C:/Dev-Cpp/lib/libgsl.a \
    C:/Dev-Cpp/lib/libgslcblas.a \
    -lm \
    C:/ASN/libs/zlib-1.2.7/libz.a \
    'C:/ASN/libs/nlopt-2.3/.libs/libnlopt_cxx.a'
INCLUDEPATH += 'C:/ASN/libs/nlopt-2.3/api' \
    'C:/Dev-Cpp/include' \
    'C:/ASN/libs/zlib-1.2.7'
QMAKE_CXXFLAGS = -O3 --fast-math
