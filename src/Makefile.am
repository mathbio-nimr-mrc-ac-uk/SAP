#===============================================================================
# Makefile.am : for automake
#===============================================================================

bin_PROGRAMS = sap

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

sap_SOURCES = \
bestrot.c bestrot.h cones.c cones.h galloc.h geom.c geom.h \
matrix.c matrix.h pdbprot.c pdbprot.h ql.c ql.h sap.c sap.h \
siva.c siva.h sort.c sort.h student.h stutest.c util.c util.h

sap_LDADD = $(INTI_LIBS)

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

