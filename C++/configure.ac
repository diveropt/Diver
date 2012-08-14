# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT([DE], [0.0.1])
AC_CONFIG_SRCDIR([src/Trial.cpp])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

LT_INIT 

# Checks for MPI.
# If --with-mpi=auto is used, try to find MPI, but use standard C compiler if it is not found.
# If --with-mpi=yes is used, try to find MPI and fail if it isn't found.
# If --with-mpi=no is used, use a standard C compiler instead.
AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi], [compile with MPI (parallelization) support. If none is found, MPI is not used. Default: auto])],,[with_mpi=auto])

AX_PROG_CXX_MPI([test x"$with_mpi" != xno],[
  ],[
  if test x"$with_mpi" = xyes; then
    AC_MSG_FAILURE([MPI compiler requested, but couldn't use MPI.])
  else
    AC_MSG_WARN([No MPI compiler found, won't use MPI.])
  fi
])


# Checks for programs.
AC_PROG_CXX_C_O
AM_INIT_AUTOMAKE

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
