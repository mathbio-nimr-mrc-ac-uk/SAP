#ifndef __WT_SORT_H__
#define __WT_SORT_H__

/*typedef struct {  int a,b; float s; char c; } Pairs;*/

typedef struct {
	int	a, b;
	float   s;
	float	c;
} Pairs;

void sorts (short *a, int n );
void sorti (int* a, int n );
void sortf ( float* a, int n );
void sort ( Pairs *pa, float *fa, int *ia, int *p, int n, int init_pointers );

#endif
