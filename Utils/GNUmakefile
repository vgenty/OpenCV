#
# This is an example GNUmakefile for my packages
#

# specific names for this package
DICT  = OpenCV_UtilsDict
SHLIB = libOpenCV_Utils.so
SOURCES = $(filter-out $(DICT).cxx, $(wildcard *.cxx))
FMWK_HEADERS = LinkDef.h $(DICT).h NDArrayConverter.h
HEADERS = $(filter-out $(FMWK_HEADERS), $(wildcard *.h))
OBJECTS = $(SOURCES:.cxx=.o)

# include options for this package
INCFLAGS  = -I.                       #Include itself
INCFLAGS += $(shell python-config --includes)
INCFLAGS += -I$(shell python -c "import numpy; print numpy.get_include()")
INCFLAGS += -I/usr/local/include/

# platform-specific options
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

LDFLAGS += -L$(shell python-config --prefix)/lib/
LDFLAGS += $(shell python-config --ldflags)
LDFLAGS += -L/usr/local/lib/ -lopencv_core -lopencv_highgui -lopencv_imgcodecs


# call kernel specific compiler setup
include $(LARLITE_BASEDIR)/Makefile/Makefile.${OSNAME}

# call the common GNUmakefile
include $(LARLITE_BASEDIR)/Makefile/GNUmakefile.CORE
