#ifndef __STUDENT_H__
#define __STUDENT_H__

/* ==== HEADER student.h ==== */

/* t-test and F-test: originally from Numerical Recipes
 * Converted into ANSI C from a Turbo Pascal 5.0 unit
 * originally written in Hungary, 7-Nov-91
 */

/* ANSI C, IRIX 5.2, 5. Oct. 1994. Andris */

/* ---- STANDARD HEADERS ---- */


/* ---- PROTOTYPES ---- */

void stutest (double Avg1, double Avg2, double Var1, double Var2, int n1, int n2);

/* t_test: Student's t-test for two distributions that have the same
 * "true" variance. Ave1, 2 are the averages, Var1, 2 are the variances.
 * (These should pass an F-test for the Var1==Var2 hypothesis.)
 * n1, n2 are the number of data points. 
 * Return values: returns the probability (significance) level 
 * and also returns the t-stat value in *T if T!=NULL. Returns -1.0 if
 * something was fishy plus prints a warning on stderr.
 */
double t_test(double Ave1, double Ave2, double Var1, double Var2, 
	int n1, int n2, double *T);

/* tu_test: Student's t-test for two distributions with different
 * variances. Ave1, 2 are the averages, Var1, 2 are the variances.
 * (These should pass an F-test for the Var1!=Var2 hypothesis.)
 * n1, n2 are the number of data points. 
 * Return values: returns the probability (significance) level 
 * and also returns the t-stat value in *T if T!=NULL. Returns -1.0 if
 * something was fishy plus prints a warning on stderr.
 */
double tu_test(double Ave1, double Ave2, double Var1, double Var2, 
	int n1, int n2, double *T);

/* f_test: Fischer's F-test to decide whether Var1, 2 are different
 * variances. n1, n2 are the no. of data points.
 * Return values: returns the F-statistics probability (significance) 
 * level. Also returns the F-value in *Fval if Fval!=NULL.
 * Returns -1.0 if something silly happens.
 */
double f_test(double Var1, double Var2, int n1, int n2, double *Fval);

/* ==== END OF HEADER student.h ==== */

#endif
