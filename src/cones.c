/*
 * Conic accessibilities (beta=sidechain centroid) on the PDB.
 * Monomeric proteins only,  output in xmgr format.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/* ---- INCLUDE FILES ---- */

#include "matrix.h"
#include "galloc.h"
#include "cones.h"

/* ANSI C, IRIX 4.0.5, 23. May 1994. Andris */

/* ---- STANDARD HEADERS ---- */

/*static int get_names(const char *Names, char Dsspnm[], char Pdbnm[]);*/

static int get_cacentr(const Chain_ * Chain, Atom_ ** Calpha,
		       Atom_ ** Centr);
static void make_distmats(const Atom_ Calpha[], const Atom_ Centr[],
			  int Rno, Trimat_ Dista, Sqmat_ Distab,
			  Trimat_ Distb);
static void betacone_shield(const Trimat_ Dista, const Sqmat_ Distab,
			    const Trimat_ Distb, float Relshb[], int Rno);



void cones(Pdbentry_ * Entry)
{
    float *Relsh = NULL;
    Trimat_ Dista, Distb;
    Sqmat_ Distab;
    Atom_ *Calpha = NULL, *Centr = NULL;
    int i, Len;

    get_cacentr(Entry->Chains, &Calpha, &Centr);
    Len = Entry->Chains->Aano;
    Dista = alloc_trimat(Len);
    assert(Dista);
    Distb = alloc_trimat(Len);
    assert(Distb);
    Distab = alloc_sqmat(Len);
    assert(Distab);
    make_distmats(Calpha, Centr, Len, Dista, Distab, Distb);
    free(Calpha);
    free(Centr);
    Calpha = Centr = NULL;
    Relsh = (float *) calloc(Len, sizeof(float));
    assert(Relsh);
    betacone_shield(Dista, Distab, Distb, Relsh, Len);
    free_matrix(Dista, Len);
    free_matrix(Distb, Len);
    free_matrix(Distab, Len);
    for (i = 0; i < Len; i++)
	Entry->Chains->Atoms[i].Bfact = Relsh[i];
    /* cleanup */
    free(Relsh);
    Relsh = NULL;
}

/* ==== FUNCTIONS ==== */

/* get_names: constructs the full DSSP and PDB pathnames out of a file
 * (Names) that contains 4-char PDB names (like 4INS) a <newline>
 * and nothing else in every line. When get_names() is called
 * for the first time, then Names is opened and the first code is read
 * and two filenames, decorated with appropriate pathnames and extensions for
 * the DSSP and PDB databases are copied into Dsspnm[] and Pdbnm[]
 * (which should be long enough to hold these) and 1 is
 * returned. The following names are returned on subsequent calls.
 * When Names is exhausted, it is closed and 0 is returned.
 * 0 is also returned in case of an I/O error.
 */
#define NAMELEN 4
/*static int get_names(const char *Names, char Dsspnm[], char Pdbnm[])
{
    static FILE *Nf = NULL;
    char Name[NAMELEN + 1];
    static const char *Dsspath = "/nonhom/dssp/", *Pdbpath = "/nonhom/",
	*Ext = ".dssp";

    if (Nf == NULL)
	Nf = fopen(Names, "r");
    if (Nf == NULL)
	return (0);
    if (1 == fscanf(Nf, "%4s\n", Name)) {
	Name[NAMELEN] = '\0';
	strcpy(Dsspnm, Dsspath);
	strcpy(Pdbnm, Pdbpath);
	strcat(Dsspnm, Name);
	strcat(Pdbnm, Name);
	strcat(Dsspnm, Ext)
	return (1);
    } else {
	fclose(Nf);
	Nf = NULL;
	return (0);
    }
}
*/
#undef NAMELEN
/* END of get_names */

/* ---- AB CONIC SHIELDEDNESS ---- */

/* get_cacentr: extracts the C-alpha coordinates from a PDB chain
 * pointed to by Chain and copies them into an array allocated
 * within and pointed to by Calpha. Also calculates the centroids
 * of the side chains and puts them into an array (also alloc-d within)
 * pointed to by Centr. If no side chain is found (e.g. Gly) then
 * the corresponding C-alpha coordinate is used as side chain centroid.
 * Both arrays are Chain->Aano long (this value is returned).
 * Returns 0 if Chain is not a protein (in this case the arrays
 * will be NULL). Revised 23-May-94 to take care of Alt-s.
 */
static int get_cacentr(const Chain_ * Chain, Atom_ ** Calpha,
		       Atom_ ** Centr)
{
    Atom_ *Cur, *Ca, *Sctr;
    int i, j, Scno = 0, Resno = -9999;
    char Rid = '\0';

    *Calpha = *Centr = NULL;
    if (Chain->Aano <= 0 || Chain->Type == 'X')
	return (0);

    Ca = (Atom_ *) calloc(Chain->Aano, sizeof(Atom_));
    assert(Ca);
    Sctr = (Atom_ *) calloc(Chain->Aano, sizeof(Atom_));
    assert(Sctr);

    /* trundle along the chain: the C-alpha always comes
     * before the corresponding side chain (we hope)
     * Check for Alt conformations: a new CA is accepted
     * only if its Resno or Rid field _changes_ compared
     * to those of the previous one. The alternative sidechain
     * conformations are lumped into one avg chain.
     */
    for (Cur = Chain->Atoms, i = j = 0; i < Chain->Atomno; i++, Cur++) {
	/* skip main-chain atoms (except C-alphas), terminal
	 * carboxyl-O-s and everybody else not C, N, O, S
	 */
	if (!strcmp(Cur->Id, "N") || !strcmp(Cur->Id, "C") ||
	    !strcmp(Cur->Id, "O") || !strcmp(Cur->Id, "OXT") ||
	    NULL == strchr("CNOS", Cur->Id[0]))
	    continue;

	/* finish up prev sidechain if any, copy the new C-alpha, 
	 * reset the centroid coords, store the AA type char
	 */
	if (!strcmp(Cur->Id, "CA")) {
	    if (Cur->Resno != Resno || Cur->Rid != Rid) {
		Resno = Cur->Resno;
		Rid = Cur->Rid;
		Ca[j] = *Cur;
		Sctr[j].X = Sctr[j].Y = Sctr[j].Z = 0.0;
		if (Scno) {
		    Sctr[j - 1].X /= Scno;
		    Sctr[j - 1].Y /= Scno;
		    Sctr[j - 1].Z /= Scno;
		    Scno = 0;
		} else if (j)
		    Sctr[j - 1] = Ca[j - 1];	/* copy C-alpha */
		j++;
	    }
	    continue;
	}

	/* side chain atoms: average in the centroid */
	Sctr[j - 1].X += Cur->X;
	Sctr[j - 1].Y += Cur->Y;
	Sctr[j - 1].Z += Cur->Z;
	Scno++;
    }				/* for Cur */

    /* finish up last */
    if (Scno) {
	Sctr[j - 1].X /= Scno;
	Sctr[j - 1].Y /= Scno;
	Sctr[j - 1].Z /= Scno;
    } else
	Sctr[j - 1] = Ca[j - 1];	/* copy C-alpha */

    *Calpha = Ca;
    *Centr = Sctr;
    return (Chain->Aano);
}

/* END of get_cacentr */

/* make_distmats: calculates the (squared) C-alpha:C-alpha,
 * C-alpha:sidechain-centroid and sidechain-centroid:sidechain-centroid
 * distances and puts them into Dista, Distab and Distb, 
 * respectively. The coordinates are in Calpha[] and Centr[], 
 * both arrays are Rno long. The matrices are assumed to be
 * allocated accordingly.
 */
static void make_distmats(const Atom_ Calpha[], const Atom_ Centr[],
			  int Rno, Trimat_ Dista, Sqmat_ Distab,
			  Trimat_ Distb)
{
    register int i, j;
    double D;

    for (i = 0; i < Rno; i++) {
	/* Ca-Ca and Centr-Centr distances */
	for (j = 0; j < i; j++) {
	    D = atom_dist(Calpha + i, Calpha + j);
	    Dista[i][j] = D * D;
	    D = atom_dist(Centr + i, Centr + j);
	    Distb[i][j] = D * D;
	}

	/* Ca[i]:Centr[j] distances */
	for (j = 0; j < Rno; j++) {
	    D = atom_dist(Calpha + i, Centr + j);
	    Distab[i][j] = D * D;
	}
    }
}

/* END of make_distmats */

/* betacone_shield: calculates the local shieldedness values for 
 * each sidechain centroid and puts them into Relshb[]. 
 * This is the "traditional" algorithm. For the k-th
 * point, all points which are closer than NBRADIUS are
 * selected, their distances from their common local centroid
 * is calculated and the angle of the smallest cone that encompasses the
 * whole set and centred on 'k' with an axis going through
 * the centroid is determined.
 */
#define NBRADIUS (8.0)
#define NBRADIUS2 NBRADIUS*NBRADIUS
static void betacone_shield(const Trimat_ Dista, const Sqmat_ Distab,
			    const Trimat_ Distb, float Relshb[], int Rno)
{
    register int i, j, k, ci, cj, Closeno;
    double *Di0, *Dik;
    int *Close;
    register double Ang, Largang, D, Trisum, Isum, Rsh;

    /* init storage for distances and indices of "close" points,
     * The "canonical ordering"
     * here is that 0..Rno-1 contains the BETAs, and Rno..2*Rno-1
     * the ALPHAs. This way the three dist matrices can be joined
     * into one big overall distmat without any clumsy indexing
     * tricks. Re-allocated every time (->static in DRAGON)
     */
    Di0 = (double *) calloc(2 * Rno, sizeof(double));	/* centroid dist sq */
    assert(Di0);
    Dik = (double *) calloc(2 * Rno, sizeof(double));	/* i-k dist sq */
    assert(Dik);
    Close = (int *) calloc(2 * Rno, sizeof(int));	/* index lookup */
    assert(Close);

    /* scan all sidechain ("beta") points */
    for (k = 0; k < Rno; k++) {
	/* select close points: an index < Rno means the index-th beta,
	 * index>=Rno means the (index-Rno)-th alpha
	 */
	Trisum = 0.0;
	for (Closeno = 0, i = 0; i < 2 * Rno; i++) {
	    if (i < Rno)	/* beta-beta */
		D = (i > k) ? Distb[i][k] : Distb[k][i];
	    else		/* beta-alpha */
		D = Distab[i - Rno][k];
	    if (D > NBRADIUS2)
		continue;	/* too far away from k */
	    Close[Closeno] = i;	/* store index */
	    Dik[Closeno++] = D;	/* store D(i,k)^2 */
	    Trisum += D;	/* start summing for local centroid */
	}

	if (Closeno <= 1) {	/* too few points, make it very exposed */
	    Relshb[k] = -1.0;
	    continue;
	}

	/* calc the distances from the local centroid using
	 * Lagrange's Theorem. The first step is to sum all
	 * interpoint distances: the i-k distances were done
	 * in the previous cycle. Note that if i<j then
	 * Close[i]<Close[j]
	 */
	for (i = 0; i < Closeno; i++)
	    for (j = 0; j < i; j++) {
		ci = Close[i];
		cj = Close[j];
		if (ci < Rno && cj < Rno)	/* beta-beta */
		    D = (ci > cj) ? Distb[ci][cj] : Distb[cj][ci];
		else if (ci < Rno && cj >= Rno)
		    D = Distab[cj - Rno][ci];
		else if (ci >= Rno && cj < Rno)
		    D = Distab[ci - Rno][cj];
		else {		/* alpha-alpha */

		    ci -= Rno;
		    cj -= Rno;
		    D = (ci > cj) ? Dista[ci][cj] : Dista[cj][ci];
		}
		Trisum += D;
	    }
	Trisum /= (Closeno + 1) * (Closeno + 1);

	/* now get squared distances for the i-th point
	 * from the centroid. If the local dist set is
	 * non-metric enough then this dist may be negative;
	 * cheat by taking the abs value
	 */
	for (i = 0; i < Closeno; i++) {
	    Isum = Dik[i];
	    ci = Close[i];
	    for (j = 0; j < Closeno; j++) {
		cj = Close[j];
		if (ci < Rno && cj < Rno)	/* beta-beta */
		    D = (ci > cj) ? Distb[ci][cj] : Distb[cj][ci];
		else if (ci < Rno && cj >= Rno)
		    D = Distab[cj - Rno][ci];
		else if (ci >= Rno && cj < Rno)
		    D = Distab[ci - Rno][cj];
		else {		/* alpha-alpha */

		    ci -= Rno;
		    cj -= Rno;
		    D = (ci > cj) ? Dista[ci][cj] : Dista[cj][ci];
		    ci += Rno;	/* restore */
		}
		Isum += D;
	    }
	    Di0[i] = fabs(Isum / (Closeno + 1) - Trisum);
	}

	/* get the dist of the k-th point from the centroid
	 * in the same way and put into D
	 */
	D = 0.0;
	for (i = 0; i < Closeno; i++)
	    D += Dik[i];
	D = fabs(D / (Closeno + 1) - Trisum);

	/* calc angle using the Cosine Rule for each entry in
	 * Close[] and determine the maximum
	 */
	Largang = -1000.0;
	for (i = 0; i < Closeno; i++) {	/* scan all angles */
	    Ang = acos((Dik[i] + D - Di0[i]) / (2.0 * sqrt(D * Dik[i])));
	    if (Ang > Largang)
		Largang = Ang;
	}

	/* save the rel. shield of largest angle for the k-th point */
	Rsh = (Largang - M_PI_2) / M_PI_2;
	if (fabs(Rsh) > 1.0) {
	    fprintf(stderr, "? betacone_shield(): Relsh=%.2e\n", Rsh);
	    Rsh = 0.0;		/* cheat again */
	}
	Relshb[k] = Rsh;
    }				/* for k */

    /* cleanup */
    free(Di0);
    free(Dik);
    free(Close);
}

#undef NBRADIUS2
#undef NBRADIUS
/* END of betacone_shield */

/* ==== END OF PROGRAM dsspabcones.c ==== */
