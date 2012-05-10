#ifndef SAPIT_H_
#define SAPIT_H_


#define NALLOC 1000
#define NACID 30

typedef struct {
	Vec	v;
	float	d;
	Vec	cos;
} Tri;

typedef struct {
	int	c, r;
} Cel;

/*typedef struct {
	int	a, b;
	float	c;
} Pairs;*/

typedef struct {
	char	*res;
	float	*acc;
	Vec	*ca, *cb;
	int	len;
} Seq;

void stats (int half, float **data, int n);
int selsort (const void *ac, const void *bc);
float compare (int cycles, Tri **a, Tri **b, Cel **c, Cel **d, Seq *seqa, Seq *seqb, int print);
float local_rms(int m, int n, int lena, int lenb, Tri **a, Tri **b);
float recycle (int cycle, Seq *seqa, Seq *seqb, float **bias, Pairs *sel, float **sec, float **sim, Tri **a, Tri **b, int *nsel, int print);
int cellhits (Cel *lista, Cel *listb, Tri **a, Tri **b, int m, int n, int lena, int lenb);
int check_sel (int **aln, Pairs *sel, int nsel, int naln);
void trace_mat (int **mat, int lena, int lenb);
void print_mat (float scale, float **mat, int lena, int lenb);
void score_pair (float dif_wt, float **bias, Pairs *sel, float **sim, Tri **a, Tri **b, int la, int lb, int nsel);
void score_pair (float dif_wt, float **bias, Pairs *sel, float **sim, Tri **a, Tri **b, int la, int lb, int nsel);
void score (float wt, float **bias, float **sim, Tri **a, Tri **b, int m, int n, int la, int lb);
float rescore (Tri **a, Tri **b, int m, int n, int **aln, int len);
float add_path (float **sim, float **smn, int na, int nb, int m, int n);
float get_path (int **aln, float **sim, int na, int nb, int *length);
int trace (float **s, int **p, int **a, int n, int i, int j);
int protin (Pdbentry_ *prot, Seq *seq, Tri ***m, Cel ***c, float z, int flip);
void add_cb (Seq *seq);
int celsort (const void *ac, const void *bc);
void set_vect (Vec *a, Vec *b, Tri **m, Cel **c, int l);
void set_cbcb (Vec *a, Vec *b, Tri **m, int l);
void extend (Vec *res, int i, int j, int k, int new);
int copyca (Chain_ *pdb, Seq *s, int flip, float z);
void flipseq (Vec *ca, char *seq, float *acc, int n);
int getca (Vec *res, FILE *pdb);
void putpdb (Seq *seq, FILE *out, char id);
void setframe (Vec a, Vec b, Vec c, Mat *frame);
int norm (float sigcut, float *data, int n);
int normn (float sigcut, float **data, int m, int n);
void matin();
void oldmatin(char *file, int mat[NACID][NACID]);
void moment (float **mom, Vec *struc, float *weight, int natom);
float super (Tri **stra, Tri **strb, Seq *seqa, Seq *seqb, float **sim, int **aln, int len);
void superout (Seq *seqa, Seq *seqb, float *sa, float *sb, int **aln, int len);
void matplot (char *file, float **mat, int lena, int lenb, int logs);


#endif /*SAPIT_H_*/
