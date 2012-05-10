#include <stdio.h>
#include <math.h>
#include "util.h"

unsigned int pack(int i, int j, int n)
/*	references top half (+diag) of a 2D array (N*N) as a vector	*/
{
    unsigned int ip, jp;
    if (n <= 0) {
	printf("Bad dimension in IPACK\n");
	return -1;
    }
    if (i <= 0) {
	printf("Bad I in IPACK: i=%d\n", i);
	return -2;
    }
    if (j <= 0) {
	printf("Bad J in IPACK: j=%d\n", j);
	return -3;
    }
    ip = max(i, j);
    jp = min(i, j) - 1;
    return (unsigned int) (ip - jp + n * jp - (jp * jp - jp) / 2);
}

void unpack(int *ip, int *jp, int n, unsigned int id)
/*	references top half (+diag) of a 2D array (N*N) from a vector	*/
{
    unsigned int i, j;
    double b = n + 0.5, a = b * b - 2 * id;
    if (n > 65535)
	printf("Error in UNPACK: n=%d is too big\n", n);
    if (a <= 0.0)
	printf("Error in UNPACK: n=%d, id=%d\n", n, id);
    j = (unsigned int) (b - sqrt(a) - 0.00001);
    i = j + id - n * j + (j * j - j) / 2;
    j++;
    if (j <= 0 || j > n)
	printf("Bad J in UNPACK: j=%d (n=%d id=%d)\n", j, n, id);
    if (i <= 0 || i > n)
	printf("Bad I in UNPACK: i=%d (n=%d id=%d)\n", i, n, id);
    *ip = min(i, j);
    *jp = max(i, j);
}

int read_line(FILE * file, char *string)
{
    char c;
    int i = 0;
    *string = 0;
    while ((c = getc(file))) {
	/* printf("%d >%c<\n", c,c); */
	if (feof(file))
	    return -i - 1;
	if (c == '\n')
	    return i;
	string[i] = c;
	string[++i] = 0;
    }
    printf("Error in file reading!\n");
    return (-1);
}

int next_line(FILE * file)
{
    char c;
    while ((c = getc(file))) {
	if (feof(file))
	    return 0;
	if (c == '\n')
	    return 1;
    }
    printf("Error in file reading!\n");
    return (0);
}

int min(int i, int j)
{
    if (i < j)
	return i;
    else
	return j;
}

int max(int i, int j)
{
    if (i > j)
	return i;
    else
	return j;
}

/* float fmax(float i,float j) { if(i>j) return i; else return j; } */
