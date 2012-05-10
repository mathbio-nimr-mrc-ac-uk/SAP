#ifndef __MATRIX_H__
#define __MATRIX_H__

/* ==== HEADER matrix.h ==== */

/* Header for square and lower triangle matrices: the latter
 are stored economically. */

/* ANSI C, IRIX 5.2, 5. Aug. 1994. Andris Aszodi */

/* ---- HEADERS ---- */


/* ---- GLOBAL TYPES ---- */

/* WARNING: The routines do not check whether the matrix is triangular
 or square. The typedef's below might help during compilation. */

typedef double **Matrix_ ;	/* general type */
typedef double **Trimat_ ;	/* triangular */
typedef double **Sqmat_ ;	/* square */

/* ---- PROTOTYPES ---- */

Trimat_ alloc_trimat(int Size);
void free_matrix(double **Mat, int Size);
void list_trimat(Trimat_ Mat, int Size, int Linewidth,
	int Width, int Prec);
void list_matptr(double **Mat, int Rno);

Sqmat_ alloc_sqmat(int Size);
void list_sqmat(Sqmat_ Mat, int Size, int Linewidth,
	int Width, int Prec);

int lu_decomp(Sqmat_ A, int n, int **Perm);
double lu_det(const Sqmat_ Lu, int Psign, int n);
void lu_solve(const Sqmat_ A, const int Perm[], double b[], int n);

/* ==== END OF HEADER matrix.h ==== */

#endif		/* __MATRIX_H__ */
