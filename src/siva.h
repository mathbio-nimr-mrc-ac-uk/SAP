#ifndef __SIVA_H__
#define __SIVA_H__

/* ==== HEADER siva.h ==== */

/* Singular value decomposition: homebrew version */

/* ANSI C, IRIX 4.0.5, 21. July 1994. Andris Aszodi */



/* ---- DEFINITIONS ---- */

#ifdef DBL_EPSILON	/* precision */
#define SIVA_EPSILON (10.0*DBL_EPSILON)
#else
#define SIVA_EPSILON (1.0e-10)
#endif

#ifdef FLT_MIN	    /* smallest acceptable eigenvalue */
#define SIVA_MINVAL (100.0*FLT_MIN)
#else
#define SIVA_MINVAL (1.0e-37)	/* required by ANSI */
#endif

/* ---- PROTOTYPES ---- */

int siva_setup(int Row, int Col, double ***U, double **W, double ***V);
int siva_decomp(const double **A, int Row, int Col, 
	double **U, double *W, double **V);
int rank_cond(double W[], int N, double Eps, double *Cond);
void siva_solve(const double **U, const double W[], const double **V, 
	int Row, int Col, const double B[], double X[]);
void free_siva(int Row, int Col, double **U, double *W, double **V);

/* ==== END OF HEADER siva.h ==== */

#endif
