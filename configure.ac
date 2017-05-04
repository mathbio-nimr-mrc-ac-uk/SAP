#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_INIT(sap, 1.1.3, jens.kleinjung@crick.ac.uk)

AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([src/sap.c])
AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Checks for programs. 
AC_PROG_CC_C_O 
AC_PROG_INSTALL 
 
# Checks for libraries. 
AC_CHECK_LIB([m], [main], [], 
                AC_MSG_ERROR([Could not find required C math library.])) 
 
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
AC_CONFIG_HEADERS([src/config.h])
AC_OUTPUT

