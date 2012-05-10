#ifndef __GALLOC_H__
#define __GALLOC_H__

/* ==== HEADER galloc.h ==== */

/* General memory allocation macros that make life easier. */

/* ANSI C, IRIX 4.0.5, 28. Jan. 1994. Andris */

/* ---- STANDARD HEADERS ---- */

#include <malloc.h>

/* ---- MACROS ---- */

#define G_MALLOC(PTR, TYPE) (PTR)=(TYPE*) malloc(sizeof(TYPE))
#define G_CALLOC(PTR, TYPE, SIZE) (PTR)=(TYPE*) calloc((SIZE), sizeof(TYPE))
#define G_REALLOC(PTR, TYPE, SIZE) \
    (PTR)=(TYPE*) realloc((PTR), (SIZE)*sizeof(TYPE))

#define GFREE(PTR) if ((PTR)!=NULL) { free(PTR); (PTR)=NULL; }

/* ==== END OF HEADER galloc.h ==== */
#endif
