/* ==== FUNCTIONS siva.c ==== */

/* Singular value decomposition: homebrew version */

/* ANSI C, IRIX 4.0.5, 21. July 1994. Andris Aszodi */

/* ----	HEADER AND INCLUDE FILES ---- */
/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <assert.h>

#include "matrix.h"
#include "ql.h"
#include "siva.h"

/* ==== FUNCTIONS ==== */

/* ---- SINGULAR VALUE DECOMPOSITION AND BACK-SUBSTITUTION ---- */

/* siva_setup: allocates the matrices etc. for storing the
 * decomposed form. The matrix to be decomposed is (Row x Col), 
 * where Row>=Col. If Row<Col, then the rows in U will be padded
 * to get a Col x Col matrix. The actual Row is returned.
 * U is (Row x Col), W is Col long, V is (Col x Col).
 */
int siva_setup(int Row, int Col, double ***U, double **W, double ***V)
{
    double **u, *w, **v;	/* local matrices and vectors */
    int i;

    if (Row <= 0 || Col <= 0) {	/* para */
	fprintf(stderr, "? siva_setup(): Bad dimensions (%d x %d)\n",
		Row, Col);
	return (Row);
    }

    /* If Row<Col, then u is padded with Col-Row zeroed rows 
     * to make it a Col x Col matrix 
     * and a warning is printed. 
     */
    if (Row < Col) {
	fprintf(stderr, "? siva_setup(): %d x %d matrix, rows padded\n",
		Row, Col);
	Row = Col;
    }
    u = (double **) calloc(Row, sizeof(double *));
    assert(u);
    for (i = 0; i < Row; i++) {
	u[i] = (double *) calloc(Col, sizeof(double));
	assert(u[i]);
    }

    w = (double *) calloc(Col, sizeof(double));
    assert(w);
    v = (double **) calloc(Col, sizeof(double *));
    assert(v);
    for (i = 0; i < Col; i++) {
	v[i] = (double *) calloc(Col, sizeof(double));
	assert(v[i]);
    }

    /* return */
    *U = u;
    *W = w;
    *V = v;
    return (Row);		/* actual no. of rows, may be padded */
}

/* END of siva_setup */

/* siva_decomp: singular value decomposition routine 
 * based on a Hungarian linear algebra book by Pa'l Ro'zsa.
 * It is assumed that
 * siva_setup() has been called prior to siva_decomp(). The "padding"
 * is taken care of automatically so Row<Col entries are accepted.
 * (DO NOT use the Row value returned by siva_setup()! )
 * SVD is carried out on A (which will be preserved).
 * The result of the SVD A=UWV' 
 * is returned in the arrays allocated by siva_setup():
 * arrays: U is a Row x Col (or Col x Col) matrix, W contains the
 * singular values in a vector and V is a Col x Col matrix.
 * Return value: 1 if the built-in iteration limit in eigen_ql()
 * is exceeded, 0 if OK.
 */
int siva_decomp(const double **A, int Row, int Col,
		double **U, double *W, double **V)
{
    Trimat_ Ata = NULL;
    register double Temp;
    register int i, j, k;
    int Err = 0;

    /* construct the A'A matrix lower triangle in Ata */
    Ata = alloc_trimat(Col);
    assert(Ata);
    for (i = 0; i < Col; i++)
	for (j = 0; j <= i; j++) {
	    Temp = 0.0;
	    for (k = 0; k < Row; k++)
		Temp += A[k][i] * A[k][j];
	    Ata[i][j] = Temp;
	}

    /* get the eigenvalues and eigenvectors of Ata:
     * W holds the eigenvalues sorted in decreasing
     * order, V holds the corresponding eigenvectors
     * as ROWS (also sorted)
     */
    Err = eigen_ql(Ata, Col, W, V);
    free_matrix(Ata, Col);
    if (Err)
	return (1);

    /* transpose V in place and sqrt the singular value vector W */
    for (i = 0; i < Col; i++) {
	Temp = W[i];
	W[i] = (Temp < SIVA_MINVAL) ? 0.0 : sqrt(Temp);
	for (j = 0; j < i; j++) {
	    Temp = V[i][j];
	    V[i][j] = V[j][i];
	    V[j][i] = Temp;
	}
    }

    /* get the matrix U: Ro'zsa says that A*v(j)=W[j]*u(j),
     * where u(j) and v(j) are the j-th columns of U and V, 
     * respectively. If W[j] is too small, then the corresponding
     * u(j) will be set to 0 (cf. sqrt(W) above)
     */
    for (j = 0; j < Col; j++) {
	if (W[j] == 0.0) {	/* set U[][j] to zero */
	    for (i = 0; i < Row; i++)
		U[i][j] = 0.0;
	    continue;
	}

	for (i = 0; i < Row; i++) {
	    Temp = 0.0;
	    for (k = 0; k < Col; k++)
		Temp += A[i][k] * V[k][j];
	    U[i][j] = Temp / W[j];
	}

	/* If Row<Col, the rest of the U rows should be zero.
	 * siva_setup() does this but as a precaution we
	 * zero these rows here anyway
	 */
	if (Row < Col)
	    for (i = Row; i < Col; i++)
		memset(U[i], 0, Col * sizeof(double));

    }
    return (0);
}

/* END of siva_decomp */

/* rank_cond: checks the N singular values W[] of a matrix 
 * after SVD. If Cond!=NULL, then the condition number 
 * (ratio of the largest and smallest singular value) is
 * calculated. The singular values which are smaller than
 * Eps times the largest are set to 0.0.
 * Return value: the rank of the matrix.
 */
int rank_cond(double W[], int N, double Eps, double *Cond)
{
    double Wmax = -HUGE_VAL, Wmin = HUGE_VAL;
    int i, Rank;

    /* get the largest and smallest singular value */
    for (i = 0; i < N; i++) {
	if (W[i] > Wmax)
	    Wmax = W[i];
	if (W[i] < Wmin)
	    Wmin = W[i];
    }

    /* calc the condition number: set to HUGE_VAL if Wmin==0.0 */
    if (Cond != NULL)
	*Cond = (Wmin == 0.0) ? HUGE_VAL : Wmax / Wmin;

    /* set all singular values which are smaller than Eps*Wmax
     * to zero: this is the conditioning. Calc the rank
     */
    Wmax *= fabs(Eps);
    Rank = N;
    for (i = 0; i < N; i++)
	if (W[i] < Wmax) {
	    W[i] = 0.0;
	    Rank--;
	}
    return (Rank);
}

/* END of rank_cond */

/* siva_solve: back-substitution routine for solving linear equations
 * AX=B. A should be SV-decomposed into U, W and V' by siva_comp()
 * and the weight vector should be "conditioned" (small entries
 * zeroed) by rank_cond() prior to the call to this routine.
 * U, W and V' are the decomposed and conditioned bits, Row and Col
 * are the row and column numbers of A (Row may have been "padded"
 * by siva_comp()!), B[] is the Row long right-hand-side vector.
 * (Padding may be necessary here, too!)
 * The result is returned in X[] (Col long).
 */
void siva_solve(const double **U, const double W[], const double **V,
		int Row, int Col, const double B[], double X[])
{
    register double *Tmp = NULL;
    register double Sum;
    register int i, j;

    /* get WU'*B first */
    Tmp = (double *) calloc(Col, sizeof(double));	/* calloc() zeroing is essential here */
    assert(Tmp);
    for (j = 0; j < Col; j++)
	if (W[j] != 0.0) {	/* skip zeroed */
	    Sum = 0.0;
	    for (i = 0; i < Row; i++)
		Sum += U[i][j] * B[i];
	    Tmp[j] = Sum / W[j];
	}

    /* multiply Tmp by V to get solution vector */
    for (i = 0; i < Col; i++) {
	Sum = 0.0;
	for (j = 0; j < Col; j++)
	    Sum += V[i][j] * Tmp[j];
	X[i] = Sum;
    }
    free(Tmp);
}

/* END of siva_solve */

/* free_siva: cleans up the space allocated to the 3 SVD arrays. */
void free_siva(int Row, int Col, double **U, double *W, double **V)
{
    int i;

    /* free up U: if Row<Col, then it is Col x Col */
    if (Row < Col)
	Row = Col;
    for (i = 0; i < Row; i++)
	free(U[i]);
    free(U);
    for (i = 0; i < Col; i++)
	free(V[i]);
    free(V);
    free(W);
}

/* END of free_siva */

/* ==== END OF FUNCTIONS siva.c ==== */
