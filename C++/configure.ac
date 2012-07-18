# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([DE], [0.0.1], [antje@fysik.su.se])
AC_CONFIG_SRCDIR([src/Trial.cpp])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

LT_INIT 
LT_LANG([C++])

# Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_C_O

# If this is Darwin make libtool replace shared object suffix ".dylib" with ".so"
AM_CONDITIONAL([FORCE_SO], [test "x$host_vendor" = xapple])

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_MEMCMP
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow sqrt strtol])

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])

AC_OUTPUT
