/* ==== FUNCTIONS ql.c ==== */

/* Eigenvalues and eigenvectors of real symmetric matrices by
 Housholder tridiagonalisation and QL-transformation.
 Adapted from Numerical Recipes. The "core" routines use 
 the 1..N indexing convention but the "shell" is ordinary 
 0..N-1 C-style. Standalone version is eigenql.c */

/* ANSI C, Iris Indigo IRIX 4.0.5, 20. Nov. 1992. Andris Aszodi */

/* ---- HEADER FILES ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "ql.h"

/* ---- TYPES ---- */

/* type of function returning int: demanded by qsort() */
typedef int (*Compfnc_) (const void *, const void *);

/* type for eigenvalue sorting */
typedef struct {
    double Eig;			/* an eigenvalue along with its */
    int Idx;			/* original position */
} Eigidx_;

/* ---- PROTOTYPES ---- */

/* eigen_ql(): prototype in "ql.h" */

static int eig_cmp(const Eigidx_ * E1, const Eigidx_ * E2);
static void tred2(double **a, int n, double *d, double *e);
static int tqli(double *d, double *e, int n, double **z, int Itno);

/* ---- FUNCTIONS ---- */

/* ---- SHELL ---- */

/* eigen_ql: a 'shell' function driving the Housholder and QL routines.
 Takes a lower triangular matrix Mat as input (with size Size), and
 produces the eigenvalues in Eval and the eigenvectors in Evec.
 Eval and Evec are assumed to be allocated outside with the
 correct sizes! Index shifts are performed
 to hack around the [1..N] convention of Numerical Recipes.
 Return value: 0 if OK, 1 if iteration limit was exceeded in tqli(). */

extern int eigen_ql(Trimat_ Mat, int Size, double *Eval, Matrix_ Evec)
{
    const int ITERNO = 30;

    int i, j, k, Err;
    Sqmat_ Qmat;
    double *Diag2;
    Eigidx_ *Evs;

    /* allocate full square matrix: input and raw eigenvectors */
    Qmat = alloc_sqmat(Size);
    assert(Qmat);

    /* copy values from lower triangle Mat */
    for (i = 0; i < Size; i++) {
	Qmat[i][i] = Mat[i][i];
	for (j = 0; j < i; j++)
	    Qmat[i][j] = Qmat[j][i] = Mat[i][j];
    }

    /* shift addresses so that [1..N] indexing be valid */
    for (i = 0; i < Size; i++)
	Qmat[i]--;
    Qmat--;

    /* create array for 2nd diagonal of tridiagonal matrix */
    Diag2 = (double *) calloc(Size, sizeof(double));
    assert(Diag2);

    /* shift their address as well for [1..N] indexing */
    Eval--;
    Diag2--;

    /* perform Housholder tridiagonalisation */
    tred2(Qmat, Size, Eval, Diag2);

    /* apply shifted QL-transforms: get eigenvalues in Diag,
       eigenvectors in Qmat */
    Err = tqli(Eval, Diag2, Size, Qmat, ITERNO);
    if (Err)
	fprintf(stderr, "Iteration limit (%d) is exceeded\n", ITERNO);

    /* shift back addresses for C-style indexing */
    Eval++;
    Diag2++;
    free(Diag2);		/* not even needed */
    for (i = 1; i <= Size; i++)
	Qmat[i]++;
    Qmat++;			/* raw eigenvectors in columns */

    /* sort eigenvalues in decreasing order */
    Evs = (Eigidx_ *) calloc(Size, sizeof(Eigidx_));
    assert(Evs);
    for (i = 0; i < Size; i++) {
	/* copy eigenvalues (rounded to 0.0 if necessary) and index */
	Evs[i].Eig = RND0(Eval[i]);
	Evs[i].Idx = i;
    }
    /* Quicksort from stdlib */
    qsort(Evs, Size, sizeof(Eigidx_), (Compfnc_) eig_cmp);
    /* permute eigenvectors according to their eigenvalues:
       and list them as rows */
    for (i = 0; i < Size; i++) {
	Eval[i] = Evs[i].Eig;	/* copy back eigenvalues */
	k = Evs[i].Idx;
	for (j = 0; j < Size; j++)
	    Evec[i][j] = Qmat[j][k];	/* copy eigenvectors */
    }

    /* output and cleanup */
    free(Evs);
    free_matrix(Qmat, Size);
    return (Err);
}

/* END of eigen_ql */

/* ---- AUXILIARY ROUTINES TO SHELL ---- */

/* eig_cmp: compares two Eigidx structs for qsort(). */

static int eig_cmp(const Eigidx_ * E1, const Eigidx_ * E2)
{
    return ((E2->Eig > E1->Eig) ? 1 : (E2->Eig < E1->Eig) ? -1 : 0);
}

/* END of eig_cmp */

/* ---- EIGENROUTINE CORE ---- */

/* tred2: Housholder tridiagonalisation. a is a real, symmetric matrix
 (size n*n). The main diagonal of the tridiag. output is returned in d,
 the second diagonal in e with e[1]==0.0. On return, a contains the
 transformation matrix "Q". [1..N] indexing convention is used 
 throughout!
 Algorithm and implementation from Numerical Recipes. A.A. has after-
 edited the function head to look like ANSI, and inserted 'register'
 local vars. */

static void tred2(double **a, int n, double *d, double *e)
{
    register int l, k, j, i;
    register double scale, hh, h, g, f;

    for (i = n; i >= 2; i--) {
	l = i - 1;
	h = scale = 0.0;
	if (l > 1) {
	    for (k = 1; k <= l; k++)
		scale += fabs(a[i][k]);
	    if (scale < EPSILON)
		e[i] = a[i][l];
	    else {
		for (k = 1; k <= l; k++) {
		    a[i][k] /= scale;
		    h += a[i][k] * a[i][k];
		}
		f = a[i][l];
		g = (RND0(f) > 0.0) ? -sqrt(h) : sqrt(h);
		e[i] = scale * g;
		h -= f * g;
		a[i][l] = f - g;
		f = 0.0;
		for (j = 1; j <= l; j++) {
		    /* Next statement can be omitted if eigenvectors not wanted */
		    a[j][i] = a[i][j] / h;
		    g = 0.0;
		    for (k = 1; k <= j; k++)
			g += a[j][k] * a[i][k];
		    for (k = j + 1; k <= l; k++)
			g += a[k][j] * a[i][k];
		    e[j] = g / h;
		    f += e[j] * a[i][j];
		}
		hh = f / (h + h);
		for (j = 1; j <= l; j++) {
		    f = a[i][j];
		    e[j] = g = e[j] - hh * f;
		    for (k = 1; k <= j; k++)
			a[j][k] -= (f * e[k] + g * a[i][k]);
		}
	    }
	} else
	    e[i] = a[i][l];
	d[i] = h;
    }
    /* Next statement can be omitted if eigenvectors not wanted */
    d[1] = 0.0;
    e[1] = 0.0;
    /* Contents of this loop can be omitted if eigenvectors not
       wanted except for statement d[i]=a[i][i]; */
    for (i = 1; i <= n; i++) {
	l = i - 1;
	if (RND0(d[i]) != 0.0) {	/* !=0.0 added */
	    for (j = 1; j <= l; j++) {
		g = 0.0;
		for (k = 1; k <= l; k++)
		    g += a[i][k] * a[k][j];
		for (k = 1; k <= l; k++)
		    a[k][j] -= g * a[k][i];
	    }
	}
	d[i] = (fabs(a[i][i]) < EPSILON) ? 0.0 : a[i][i];
	a[i][i] = 1.0;
	for (j = 1; j <= l; j++)
	    a[j][i] = a[i][j] = 0.0;
    }
}

/* END of tred2 */

/* tqli: QL algorithm with implicit shifts on tridiagonal matrices.
 The main diagonal is in d, the second diagonal is in e, with e[1]
 ignored. Size is n. If d and e were obtained from a general symmetric
 real matrix by Housholder transformation by tred2(), then z should
 contain the transformation matrix "Q" on input; otherwise it should
 be the unit matrix. On output, d contains the eigenvalues and z the
 eigenvectors, with the k-th column corresponding to the k-th eigenvalue.
 Itno supplies the maximum allowable no. of iterations (an addition
 by A.A.). Return value: 0 if OK, 1 if iteration limit has been
 exceeded (added by A.A. to replace the nrerror() error message function
 originally used in Numerical Recipes routines).
 [1..N] indexing convention is used throughout!
 Algorithm and implementation from Numerical Recipes. A.A. has after-
 edited the function head to look like ANSI, and inserted 'register'
 local vars. */

#define SIGN(a,b) ()

static int tqli(double *d, double *e, int n, double **z, int Itno)
{
    register int m, l, iter, i, k;
    register double s, r, ra, p, g, f, c, b;
    register float dd;

    for (i = 2; i <= n; i++)
	e[i - 1] = e[i];
    e[n] = 0.0;
    for (l = 1; l <= n; l++) {
	iter = 0;
	do {
	    for (m = l; m <= n - 1; m++) {
		dd = (float) (fabs(d[m]) + fabs(d[m + 1]));
		if ((float) fabs(e[m]) + dd == dd)
		    break;
	    }
	    if (m != l) {
		if (iter++ >= Itno)	/* too many iters */
		    return (1);
		g = (d[l + 1] - d[l]) / (2.0 * e[l]);
		r = sqrt((g * g) + 1.0);
		ra = (RND0(g) < 0.0) ? -fabs(r) : fabs(r);
		g = d[m] - d[l] + e[l] / (g + ra);
		s = c = 1.0;
		p = 0.0;
		for (i = m - 1; i >= l; i--) {
		    f = s * e[i];
		    b = c * e[i];
		    if (fabs(f) >= fabs(g)) {
			c = g / f;
			r = sqrt((c * c) + 1.0);
			e[i + 1] = f * r;
			c *= (s = 1.0 / r);
		    } else {
			s = f / g;
			r = sqrt((s * s) + 1.0);
			e[i + 1] = g * r;
			s *= (c = 1.0 / r);
		    }
		    g = d[i + 1] - p;
		    r = (d[i] - g) * s + 2.0 * c * b;
		    p = s * r;
		    d[i + 1] = g + p;
		    g = c * r - b;
		    /* Next loop can be omitted if eigenvectors not wanted */
		    for (k = 1; k <= n; k++) {
			f = z[k][i + 1];
			z[k][i + 1] = s * z[k][i] + c * f;
			z[k][i] = c * z[k][i] - s * f;
		    }
		}
		d[l] = d[l] - p;
		e[l] = g;
		e[m] = 0.0;
	    }
	} while (m != l);
    }
    return (0);
}

/* END of tqli */

/* ==== END OF FUNCTIONS ql.c ==== */
