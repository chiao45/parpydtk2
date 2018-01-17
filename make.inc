CXX = mpicxx

CPP_SRC_FILES = pydtk2/DtkWrap.hpp

PYTHON_VERSION = 3
PYTHON = /usr/bin/python$(PYTHON_VERSION)
CYTHON = cython
PYX_SRC = pydtk2/dtk2mapper.pyx
PYX_OBJ = $(PYX_SRC:.pyx=.cpp)

.SUFFIXES: .pyx
