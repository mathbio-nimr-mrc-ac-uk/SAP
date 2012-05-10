#ifndef __BESTROT_H__
#define __BESTROT_H__

#include "matrix.h"

/* ==== HEADER bestrot.h ==== */

/* An implementation of the point set alignment algorithm by
 * A. D. McLachlan. Reference:
 * McLachlan, A. D. (1979): J. Mol. Biol. 128: 49-79.
 * Replaces the buggy Kabsch rotation algorithm.
 */

/* ANSI C, IRIX 5.2, 5. Aug. 1994. Andris Aszodi */


/* ---- PROTOTYPES ---- */

double supermac (const double *, double **, double **, int , double *, double *, Sqmat_);
double *center_vectors(double **X, double *Ctr, const double *W, unsigned int Vno);
double best_rot(double **X, double **Y, const double *W,
		unsigned int Vno, Sqmat_ Transform);

/* ==== END OF HEADER bestrot.h ==== */

#endif
