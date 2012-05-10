#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(sap, 1.1.1, jkleinj@nimr.mrc.ac.uk)

AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([src/sap.c])
AC_CONFIG_HEADERS([src/config.h])

AM_INIT_AUTOMAKE

AC_SUBST(INTI_CFLAGS)
AC_SUBST(INTI_LIBS)

# Checks for programs.
AC_PROG_CC
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([m], [cos])

# Checks for header files.
AC_CHECK_HEADERS([float.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 tests/Makefile])

AC_OUTPUT
