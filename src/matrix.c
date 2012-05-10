/* ==== FUNCTIONS matrix.c ==== */

/* Routine collection for square and lower triangle matrices: the latter
 are stored economically. */

/* ANSI C, IRIX 5.2, 5. Aug 1994. Andris Aszodi */

/* ---- HEADER ---- */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>		/* for precision stuff */
#include <float.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"

/* ---- DEFINITIONS ---- */

#ifdef FLT_MIN
#define LU_EPSILON (10.0*FLT_MIN)
#else
#define LU_EPSILON (1.0e-30)
#endif

/* ==== FUNCTIONS ==== */

/* ---- LOWER TRIANGLE MATRICES ---- */

/* alloc_trimat: allocates space for a triangular matrix with Size rows.
 The triangle contains the main diagonal as well. 
 Returns the pointer to the matrix or NULL if alloc failed. */

extern Trimat_ alloc_trimat(int Size)
{
    Trimat_ Mat;
    int i;

    Mat = (double **) calloc(Size, sizeof(double *));
    if (Mat != NULL)
	for (i = 0; i < Size; i++) {
	    Mat[i] = (double *) calloc(i + 1, sizeof(double));
	    if (Mat[i] == NULL) {	/* clean up already allocated rows */
		free_matrix(Mat, i);
		return (NULL);
	    }
	}
    return (Mat);
}

/* END of alloc_trimat */

/* free_matrix: frees up space occupied by Mat. The same routine for
 triangular and square matrices. */

extern void free_matrix(double **Mat, int Size)
{
    int i;

    for (i = Size - 1; i >= 0; i--)	/* dealloc in reverse order */
	free(Mat[i]);
    free(Mat);
}

/* END of free_matrix */

/* list_trimat: lists Mat to stdout with entries occupying Width chars,
 Prec digits precision. If a row takes up more than Linewidth chars,
 then the matrix is cut up nicely. */

extern void list_trimat(Trimat_ Mat, int Size, int Linewidth,
			int Width, int Prec)
{
    char Entrys[10], Cols[8], Rows[8];	/* format strings and len determ */
    int i, j, k, Sizew, Chunk, Jbeg, Items, Ulinelen;

    sprintf(Cols, "%d", Size);	/* get printed width of Size */
    Sizew = strlen(Cols);
    if (Sizew > Width)
	Width = Sizew;
    sprintf(Entrys, "%%-%d.%df ", Width, Prec);	/* make format strs */
    sprintf(Cols, "%%-%dd ", Width);
    sprintf(Rows, "%%%dd | ", Sizew);
    Items = (Linewidth - Sizew - 3) / (Width + 1);	/* columns per chunk */

    /* main cycle: print chunks of the matrix */
    for (Jbeg = 0, Chunk = (Size - 1) / Items + 1; Chunk > 0;
	 Jbeg += Items, Chunk--) {
	/* underline length */
	Ulinelen = (Chunk > 1) ? Items * (Width + 1) + Sizew + 3 :
	    (Size - Jbeg) * (Width + 1) + Sizew + 3;
	/* print head line */
	for (k = 0; k < Sizew + 3; k++)
	    putchar(' ');
	for (j = Jbeg; j < Size && j < Jbeg + Items; j++)
	    printf(Cols, j);	/* col numbers */
	putchar('\n');
	for (k = 0; k < Ulinelen; k++)
	    putchar('-');
	putchar('\n');

	/* print chunks for all rows */
	for (i = 0; i < Size; i++) {
	    printf(Rows, i);	/* row idx */
	    for (j = Jbeg; j <= i && j < Jbeg + Items; j++)
		printf(Entrys, Mat[i][j]);
	    putchar('\n');
	}

	/* chunk separator */
	putchar('\n');
	for (k = 0; k < Ulinelen; k++)
	    putchar('=');
	puts("\n");
    }
}

/* END of list_trimat */

/* list_matptr: lists the address of the matrix (i.e. the pointer Mat
 * itself) and the row pointers (the row array addresses) to 'stderr'.
 * Provided for debugging purposes.
 */
extern void list_matptr(double **Mat, int Rno)
{
    int i;
	void *p;

	p = Mat;

    fprintf(stderr, "Matrix address=%p\n", p);
    fflush(stderr);
    for (i = 0; i < Rno; i++) {
	p = Mat[i];
	fprintf(stderr, "[%d]=%p\n", i, p);
	fflush(stderr);
    }
}

/* END of list_matptr */

/* ---- SQUARE MATRICES ---- */

/* alloc_sqmat: allocates space for a square matrix (Size*Size).
 Returns the pointer to the matrix or NULL if alloc failed. */

extern Sqmat_ alloc_sqmat(int Size)
{
    Sqmat_ Mat;
    int i;

    Mat = (double **) calloc(Size, sizeof(double *));
    if (Mat != NULL)
	for (i = 0; i < Size; i++) {
	    Mat[i] = (double *) calloc(Size, sizeof(double));
	    if (Mat[i] == NULL) {	/* clean up already allocated rows */
		printf("strange\n");
		exit(1);
		free_matrix(Mat, i);
		return (NULL);
	    }
	}
    return (Mat);
}

/* END of alloc_sqmat */

/* list_sqmat: lists Mat to stdout with entries occupying Width chars,
 Prec digits precision. If a row takes up more than Linewidth chars,
 then the matrix is cut up nicely. */

extern void list_sqmat(Sqmat_ Mat, int Size, int Linewidth,
		       int Width, int Prec)
{
    char Entrys[10], Cols[8], Rows[8];	/* format strings and len determ */
    int i, j, k, Sizew, Chunk, Jbeg, Items, Ulinelen;

    sprintf(Cols, "%d", Size);	/* get printed width of Size */
    Sizew = strlen(Cols);
    if (Sizew > Width)
	Width = Sizew;
    sprintf(Entrys, "%%-%d.%df ", Width, Prec);	/* make format strs */
    sprintf(Cols, "%%-%dd ", Width);
    sprintf(Rows, "%%%dd | ", Sizew);
    Items = (Linewidth - Sizew - 3) / (Width + 1);	/* columns per chunk */

    /* main cycle: print chunks of the matrix */
    for (Jbeg = 0, Chunk = (Size - 1) / Items + 1; Chunk > 0;
	 Jbeg += Items, Chunk--) {
	/* underline length */
	Ulinelen = (Chunk > 1) ? Items * (Width + 1) + Sizew + 3 :
	    (Size - Jbeg) * (Width + 1) + Sizew + 3;
	/* print head line */
	for (k = 0; k < Sizew + 3; k++)
	    putchar(' ');
	for (j = Jbeg; j < Size && j < Jbeg + Items; j++)
	    printf(Cols, j);	/* col numbers */
	putchar('\n');
	for (k = 0; k < Ulinelen; k++)
	    putchar('-');
	putchar('\n');

	/* print chunks for all rows */
	for (i = 0; i < Size; i++) {
	    printf(Rows, i);	/* row idx */
	    for (j = Jbeg; j < Size && j < Jbeg + Items; j++)
		printf(Entrys, Mat[i][j]);
	    putchar('\n');
	}

	/* chunk separator */
	putchar('\n');
	for (k = 0; k < Ulinelen; k++)
	    putchar('=');
	puts("\n");
    }
}

/* END of list_sqmat */

/* ---- LU-DECOMPOSITION ---- */

/* lu_decomp: performs an LU-decomposition in place on the n*n matrix
 * A. Based on partial pivoting: the row permutations are done in Perm[]
 * (allocated within!) and will be used by lu_solve().
 * If **Perm==NULL, then the permutation vector will be used
 * internally and will be freed before return. This option is
 * used when only the determinant is calculated from the LU-decomposition.
 * Return value: the sign of the determinant of the permutation
 * matrix (+/-1) or 0 if A is singular or n<=0.
 */
int lu_decomp(Sqmat_ A, int n, int **Perm)
{
    register int i, j, k, imax = 0;
    register int Psign = 1;
    int *Idx = NULL;
    register double Large, Pivot, Tmp, Tmp2;
    double *Scal = NULL;

    /* check and array initialisation */
    if (n <= 0) {
	fprintf(stderr, "? lu_decomp(): invalid matrix size %d\n", n);
	if (Perm != NULL)
	    *Perm = NULL;
	return (0);
    }
    Idx = (int *) calloc(n, sizeof(int));	/* permutation vector */
    assert(Idx);
    Scal = (double *) calloc(n, sizeof(double));	/* implicit scaling array */
    assert(Scal);

    /* get implicit scaling: if a row contains 0-s only,
     * then the matrix is singular which will be indicated
     * by setting Psign=0. Precision is controlled by
     * the constant LU_EPSILON (see Definitions above).
     */
    for (i = 0; Psign && i < n; i++) {
	Large = 0.0;
	for (j = 0; j < n; j++)
	    if ((Tmp = fabs(A[i][j])) > Large)
		Large = Tmp;
	if (Large < LU_EPSILON) {	/* (almost) singular */
	    Psign = 0;
	    break;
	}
	Scal[i] = 1.0 / Large;
    }

    /* loop over columns */
    for (j = 0; Psign && j < n; j++) {
	for (i = 0; i < j; i++) {
	    Tmp = A[i][j];
	    for (k = 0; k < i; k++)
		Tmp -= A[i][k] * A[k][j];
	    A[i][j] = Tmp;
	}

	/* find largest pivot */
	Large = 0.0;
	for (i = j; i < n; i++) {
	    Tmp = A[i][j];
	    for (k = 0; k < j; k++)
		Tmp -= A[i][k] * A[k][j];
	    A[i][j] = Tmp;

	    if ((Tmp2 = Scal[i] * fabs(Tmp)) >= Large) {	/* best so far */
		Large = Tmp2;
		imax = i;
	    }
	}

	/* interchange rows? */
	if (j != imax) {
	    for (k = 0; k < n; k++) {	/* not too efficient copy */
		Tmp = A[imax][k];
		A[imax][k] = A[j][k];
		A[j][k] = Tmp;
	    }
	    Psign *= (-1);	/* parity change */
	    Scal[imax] = Scal[j];
	}
	Idx[j] = imax;

	/* get the pivot */
	Pivot = A[j][j];
	if (fabs(Pivot) < LU_EPSILON) {	/* singularity */
	    Psign = 0;
	    break;
	}

	/* divide by the pivot */
	if (j < n - 1)
	    for (i = j + 1; i < n; i++)
		A[i][j] /= Pivot;
    }				/* for j */

    free(Scal);
    if (Perm != NULL)
	*Perm = Idx;
    else
	free(Idx);
    return (Psign);
}

/* END of lu_decomp */

/* lu_det: calculates the determinant of the n x n LU-decomposed
 * square matrix Lu. Psign is the permutation sign returned by
 * lu_decomp().
 */
double lu_det(const Sqmat_ Lu, int Psign, int n)
{
    register int i;
    register double Det, Aii;

    if (!Psign)
	return (0.0);		/* matrix is singular */

    Det = 0.0;
    for (i = 0; i < n; i++) {
	if ((Aii = Lu[i][i]) < 0.0)
	    Psign *= (-1);	/* save sign */
	Det += log(fabs(Aii));	/* sum logarithms to avoid overflow */
    }
    Det = Psign * exp(Det);
    return (Det);
}

/* END of lu_det */

/* lu_solve: solves the linear equation A*x=b (A is n*n, x,b are n long).
 * A is supposed to have been LU-decomposed by lu_decomp() above and
 * the row permutation is stored in Perm[]. b[] is the "right-hand-side"
 * vector which contains the solution on return.
 */
void lu_solve(const Sqmat_ A, const int Perm[], double b[], int n)
{
    register int i, j, ip;
    register double Tmp;

    /* permute forward */
    for (i = 0; i < n - 1; i++)
	if ((ip = Perm[i]) != i) {
	    Tmp = b[ip];
	    b[ip] = b[i];
	    b[i] = Tmp;
	}

    /* forward substitution */
    for (i = 0; i < n; i++) {
	Tmp = b[i];
	for (j = 0; j < i; j++)
	    Tmp -= A[i][j] * b[j];
	b[i] = Tmp;
    }

    /* back substitution */
    for (i = n - 1; i >= 0; i--) {
	Tmp = b[i];
	for (j = i + 1; j < n; j++)
	    Tmp -= A[i][j] * b[j];
	b[i] = Tmp / A[i][i];
    }
}

/* END of lu_solve */

#undef LU_EPSILON

/* ==== END OF FUNCTIONS matrix.c ==== */
