/* ==== PROGRAM stutest.c ==== */

/* Tests whether the means of two distributions are different. */

/* ANSI C, IRIX 5.2, 14. Oct. 1994. Andris */

/* ---- STANDARD HEADERS ---- */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/* ---- INCLUDE FILES ---- */

#include "student.h"

/* ---- DEFINITIONS ---- */

#define SIGNIF (0.05)
#define VERY_SIGNIF (0.01)

/* ---- PROTOTYPES ---- */

/*static int get_avgsd(double *Avg1, double *Avg2, double *Sd1, double *Sd2,
		     int *n1, int *n2);*/
static void list_results(double Avg1, double Avg2, double Sd1,
			 double Sd2, double Fprob, double tprob);

/* ==== MAIN ==== */

/* The program reads its input from stdin which should
 * consist of six numbers separated by whitespace and
 * end in a newline:-
 *	Avg1 Sd1 n1 Avg2 Sd2 n2 \n
 * Avg1, 2 are averages, Sd1, 2 are S.D-s, n1, n2 are the
 * no. of data points for distributions 1 and 2,  respectively.
 * Output is listed to stdout and contains the significance
 * limits for the averages and SDs as supplied by the t- and
 * F-tests. An asterisk is printed for significances P<SIGNIF, 
 * two asterisks for P<VERY_SIGNIF.
 */

void stutest(double Avg1, double Avg2, double Var1, double Var2, int n1,
	     int n2)
{
    /*    double Sd1, Sd2, F, t, Fprob, tprob;*/
    double F, t, Fprob, tprob;

    /* test variance equality */
    Fprob = f_test(Var1, Var2, n1, n2, &F);
    if (Fprob < 0.0)
	printf("Cannot perform F-test\n");

    /* perform equal-variance or unequal-variance t-test
     * depending on the result of the variance F-test 
     */
    if (Fprob < SIGNIF)
	tprob = tu_test(Avg1, Avg2, Var1, Var2, n1, n2, &t);	/* unequal */
    else
	tprob = t_test(Avg1, Avg2, Var1, Var2, n1, n2, &t);	/* equal */

    if (tprob < 0.0)
	printf("Cannot perform t-test\n");
    /* output to stdout */
    list_results(Avg1, Avg2, Var1, Var2, Fprob, tprob);
}

/* ==== FUNCTIONS ==== */

/* get_avgsd: reads the averages and SDs from stdin, in the order
 * Avg1 Sd1 n1 Avg2 Sd2 n2, separated by whitespaces,  terminated
 * by newline. Returns these values in the corresponding parameters.
 * Return value: 0 at EOF, -1 if the line could not be parsed correctly,
 * -2 for comment lines (ignored silently), 
 * 1 if OK.
 */
/*static int get_avgsd(double *Avg1, double *Avg2, double *Sd1, double *Sd2,
		     int *n1, int *n2)
{
    char Line[132];

    if (NULL == fgets(Line, 130, stdin))
	return (0);
    if (Line[0] == '#')
	return (-2);
    if (6 > sscanf(Line, "%lf %lf %d %lf %lf %d",
		   Avg1, Sd1, n1, Avg2, Sd2, n2))
	return (-1);
    else
	return (1);
}
*/
/* END of get_avgsd */

/* list_results: list to stdout neatly. */
static void list_results(double Avg1, double Avg2, double Var1,
			 double Var2, double Fprob, double tprob)
{
    double Sd1, Sd2;
    /* averages */
    printf("Avg: %.2e %s %.2e Tprob=%.2e ",
	   Avg1, (tprob < SIGNIF) ? ((Avg1 < Avg2) ? "<" : ">") : "=",
	   Avg2, tprob);
    if (tprob < SIGNIF)
	putchar('*');
    if (tprob < VERY_SIGNIF)
	putchar('*');		/* second asterisk */
    putchar('\n');

    /* SDs */
    Sd1 = sqrt(Var1);
    Sd2 = sqrt(Var2);
    printf("StD: %.2e %s %.2e Fprob=%.2e ",
	   Sd1, (Fprob < SIGNIF) ? ((Sd1 < Sd2) ? "<" : ">") : "=", Sd2,
	   Fprob);
    if (Fprob < SIGNIF)
	putchar('*');
    if (Fprob < VERY_SIGNIF)
	putchar('*');		/* second asterisk */
    putchar('\n');

}

/* END of list_results */

/* ==== END OF PROGRAM stutest.c ==== */

/* ==== FUNCTIONS student.c ==== */

/* t-test and F-test: originally from Numerical Recipes
 * Converted into ANSI C from a Turbo Pascal 5.0 unit
 * originally written in Hungary, 7-Nov-91
 */

/* ANSI C, IRIX 5.2, 5. Oct. 1994. Andris */

/* ---- PROTOTYPES ---- */

static double incompl_beta(double a, double b, double x);
static double ln_gamma(double xx);
static double beta_cf(double a, double b, double x);

/* ==== FUNCTIONS ==== */

/* t_test: Student's t-test for two distributions that have the same
 * "true" variance. Ave1, 2 are the averages, Var1, 2 are the variances.
 * (These should pass an F-test for the Var1==Var2 hypothesis.)
 * n1, n2 are the number of data points. 
 * Return values: returns the probability (significance) level 
 * and also returns the t-stat value in *T if T!=NULL. Returns -1.0 if
 * something was fishy plus prints a warning on stderr.
 */
double t_test(double Ave1, double Ave2, double Var1, double Var2,
	      int n1, int n2, double *T)
{
    int Df;
    double Svar, Tval, Prob;

    Df = n1 + n2 - 2;		/* degrees of freedom */
    if (Df <= 0) {
	fprintf(stderr, "? t_test(): Invalid (%d) degrees of freedom\n",
		Df);
	return (-1.0);
    }
    Svar = ((n1 - 1) * Var1 + (n2 - 1) * Var2) / Df;	/* pooled variance */
    Tval = fabs(Ave1 - Ave2) / sqrt(fabs(Svar * (1.0 / n1 + 1.0 / n2)));
    Prob = incompl_beta(0.5 * Df, 0.5, (double) Df / (Df + Tval * Tval));
    if (T != NULL)
	*T = Tval;
    return (Prob);
}

/* END of t_test */

/* tu_test: Student's t-test for two distributions with different
 * variances. Ave1, 2 are the averages, Var1, 2 are the variances.
 * (These should pass an F-test for the Var1!=Var2 hypothesis.)
 * n1, n2 are the number of data points. 
 * Return values: returns the probability (significance) level 
 * and also returns the t-stat value in *T if T!=NULL. Returns -1.0 if
 * something was fishy plus prints a warning on stderr.
 */
double tu_test(double Ave1, double Ave2, double Var1, double Var2,
	       int n1, int n2, double *T)
{
    double Tval, Df, Prob;

    if (n1 <= 1 || n2 <= 1) {
	fprintf(stderr, "? tu_test(): Invalid no. of data (%d, %d)\n", n1,
		n2);
	return (0.0);
    }
    Var1 /= n1;
    Var2 /= n2;
    Tval = fabs(Ave1 - Ave2) / sqrt(fabs(Var1 + Var2));	/* t-value */
    /* degrees of freedom */

    Df = (Var1 + Var2) * (Var1 + Var2) /
	(Var1 * Var1 / (n1 - 1) + Var2 * Var2 / (n2 - 1));
    Prob = incompl_beta(0.5 * Df, 0.5, Df / (Df + Tval * Tval));
    if (T != NULL)
	*T = Tval;
    return (Prob);
}

/* END of tu_test */

/* f_test: Fischer's F-test to decide whether Var1, 2 are different
 * variances. n1, n2 are the no. of data points.
 * Return values: returns the F-statistics probability (significance) 
 * level. Also returns the F-value in *Fval if Fval!=NULL.
 * Returns -1.0 if something silly happens.
 */
double f_test(double Var1, double Var2, int n1, int n2, double *Fval)
{
    int Df1, Df2;
    double F, Prob;

    /* check */
    if (n1 <= 1 || n2 <= 1) {
	fprintf(stderr, "? f_test(): Invalid no. of data (%d, %d)\n", n1,
		n2);
	return (0.0);
    }

    /* swap if Var2<Var1 */
    if (Var1 > Var2) {
	F = Var1 / Var2;
	Df1 = n1 - 1;
	Df2 = n2 - 1;
    } else {
	F = Var2 / Var1;
	Df1 = n2 - 1;
	Df2 = n1 - 1;
    }
    Prob = incompl_beta(0.5 * Df2, 0.5 * Df1, Df2 / (Df2 + Df1 * F));
    if (Prob > 1.0)
	Prob = 2.0 - Prob;
    if (Fval != NULL)
	*Fval = F;
    return (Prob);
}

/* END of f_test */

/* ---- AUXILIARIES ---- */

/* incompl_beta: the incomplete beta function approximated
 * by a continued fraction a la Numerical Recipes.
 * If x is outside the range [0..1] then -1.0 is returned.
 */
static double incompl_beta(double a, double b, double x)
{
    double bt, Ibeta;

    if (x < 0.0 || x > 1.0) {
	fprintf(stderr, "? incompl_beta(): x=%.2e is out of range\n", x);
	return (-1.0);
    }
    if (x == 0.0 || x == 1.0)
	return (x);

    bt = ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * log(x) +
	b * log(1.0 - x);
    bt = exp(bt);
    Ibeta = (x < ((a + 1.0) / (a + b + 2.0))) ?
	bt * beta_cf(a, b, x) / a : 1.0 - bt * beta_cf(b, a, 1.0 - x) / b;
    return (Ibeta);
}

/* END of incompl_beta */

/* ln_gamma: returns the value of the logarithm of the Gamma function. */
static double ln_gamma(double xx)
{
    static const double STP = 2.5066282746310005;
    static const double Coeff[6] =
	{ 76.18009172947146, -86.50532032941677, 24.01409824083091,
	-1.231739572450155, 1.208650973866179e-3, -5.395239384953e-6
    };
    double x, Tmp, Ser;
    int j;

    x = xx - 1.0;
    Tmp = x + 5.5;
    Tmp -= (x + 0.5) * log(Tmp);
    Ser = 1.0;
    for (j = 0; j <= 5; j++) {
	x += 1.0;
	Ser += Coeff[j] / x;
    }
    return (log(STP * Ser) - Tmp);
}

/* END of ln_gamma */

/* beta_cf: continued fraction approximation of the beta function. */
static double beta_cf(double a, double b, double x)
{
    static const int ITMAX = 100;
    static const double EPS = 1.0e-15;

    double tem, qap, qam, qab, em, d;
    double bz, bpp, bp, bm, az, app, am, aold, ap;
    int m;

    am = bm = az = 1.0;
    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    bz = 1.0 - qab * x / qap;

    for (m = 1; m <= ITMAX; m++) {
	em = m;
	tem = 2 * em;
	d = em * (b - m) * x / ((qam + tem) * (a + tem));
	ap = az + d * am;
	bp = bz + d * bm;
	d = -(a + em) * (qab + em) * x / ((a + tem) * (qap + tem));
	app = ap + d * az;
	bpp = bp + d * bz;
	aold = az;
	am = ap / bpp;
	bm = bp / bpp;
	az = app / bpp;
	bz = 1.0;
	if ((fabs(az - aold)) < (EPS * abs(az)))
	    return (az);
    }
    fprintf(stderr, "? beta_cf(): Cannot converge in %d steps\n", ITMAX);
    return (az);
}

/* END of beta_cf */

/* ==== END OF FUNCTIONS student.c ==== */
