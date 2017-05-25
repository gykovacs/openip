################################################################################
#
# imageMath project
#
################################################################################

TARGET = gaborTests
QT -= core
QT -= gui
CONFIG -= qt
win32: CONFIG += console


# generic build settings
INCLUDEPATH += ../../lib/openipDS
INCLUDEPATH += ../../lib/openipIO
INCLUDEPATH += ../../lib/openipSC
INCLUDEPATH += ../../lib/openipML
INCLUDEPATH += ../../lib/openipLL


# header files
#HEADERS +=


# source files
SOURCES += src/gaborTests.cc


# shared libs
SHARED_LIBS += openipDS
SHARED_LIBS += openipIO
SHARED_LIBS += openipSC
SHARED_LIBS += openipML
SHARED_LIBS += openipLL


# shared system libs
#LIBS +=


# setup everything for a shared library
include(../../../pri/openip.app.pri)


# override default options
win32: test_gsl(require)
win32: test_libpng(require)
win32: test_jpeg(require)
win32: test_tiff(require)
unix: test_libpng(require)
unix: test_jpeg(require)
unix: test_tiff(require)
unix: test_minc2(require)

test_opencl()
QMAKE_CXXFLAGS -= -fno-rtti						# we need rtti
