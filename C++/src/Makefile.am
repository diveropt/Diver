## Makefile.am for src
## Process this file with automake to produce Makefile.in

MODULE = src

SRC_DIR = $(srcdir)
INC_DIR = $(top_srcdir)/include

INC_FILES = $(INC_DIR)/global.h \
$(INC_DIR)/Trial.h \
$(INC_DIR)/DE.h \
$(INC_DIR)/Mutation.h \
$(INC_DIR)/Crossover.h \
$(INC_DIR)/Selection.h

SRC_FILES = $(SRC_DIR)/Trial.cpp \
$(SRC_DIR)/Selection.cpp

AM_CPPFLAGS = -I$(INC_DIR)

lib_LTLIBRARIES = libDE.la

libDE_la_SOURCES = $(SRC_FILES)

libDE_la_includedir = $(includedir)

libDE_la_include_HEADERS = $(INC_FILES) 

libDE_la_LDFLAGS = -shared -export-dynamic

if FORCE_SO
libDE_la_LDFLAGS += -shrext ".so"
endif