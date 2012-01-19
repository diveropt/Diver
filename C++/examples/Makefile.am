## Makefile.am for examples
## Process this file with automake to produce Makefile.in

MODULE = examples

SRC_DIR = $(srcdir)/src
INC_DIR = $(srcdir)/include

PROG_LIBS =  $(top_srcdir)/src/libDE.la

AM_CPPFLAGS = -I$(INC_DIR) \
-I$(top_srcdir)/include

AM_LDFLAGS = -shared

bin_PROGRAMS = DE_example

DE_example_SOURCES = $(SRC_DIR)/DE_example.cpp 
DE_example_includedir = $(includedir)
DE_example_include_HEADERS = $(INC_DIR)/TestFunctions.h

DE_example_LDADD = $(PROG_LIBS)

