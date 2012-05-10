#ifndef __QL_H__
#define __QL_H__

/* ==== HEADER ql.h ==== */

/* Eigenvalues and eigenvectors of real symmetric matrices by
 Housholder tridiagonalisation and QL-transformation.
 Adapted from Numerical Recipes. The "core" routines use 
 the 1..N indexing convention but the "shell" is ordinary 
 0..N-1 C-style. Standalone version is eigenql.c */

/* ANSI C, Iris Indigo IRIX 4.0.5, 20. Nov. 1992. Andris Aszodi */

/* ---- INCLUDE FILES ---- */
#include "matrix.h"


/* ---- ROUNDOFF ERROR CHECK ---- */

#define EPSILON 1.0e-10
#define RND0(x) (x=(fabs(x)<EPSILON)? 0.0: (x))

/* ---- PROTOTYPES ---- */

extern int eigen_ql(Trimat_ Mat, int Size, double *Eval, Matrix_ Evec);

/* ==== END OF HEADER ql.h ==== */

#endif	/* __QL_H__ */
