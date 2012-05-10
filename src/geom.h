#ifndef __WT_GEOM_H__
#define __WT_GEOM_H__

#define PI 3.14159265358979323846
#define twoPI PI+PI

/*      STRUCTURES FOR 3D GEOMETRY      */

typedef struct { float x, y, z; } Vec;
typedef struct { Vec A, B, C; } Mat;




void vinit (Vec *c);
void vcopy (Vec b, Vec *c);
void vunit (Vec b, Vec *c);
void vnorm (Vec *c);
void vave (Vec a, Vec b, Vec *c);
void vsum (Vec a, Vec *c);
void vadd (Vec a, Vec b, Vec *c);
void vsub (Vec a, Vec b, Vec *c);
void vmul (Vec *c, float s);
void vdiv (Vec *c, float s);
float	vdif (Vec a, Vec b);
float	vddif (Vec a, Vec b);
float	vdad (Vec a, Vec b);
float	vddad (Vec a, Vec b);
float	vsqr (Vec a);
float	vmod (Vec a);
float	vdot (Vec a, Vec b);
void vprod (Vec a, Vec b, Vec *c);
float	vtri (Vec a, Vec b, Vec c);
float	pdotp (Vec a, Vec b, Vec c, Vec d);
float	pvol (Vec a, Vec b, Vec c, Vec d);
float	phand (Vec a, Vec b, Vec c, Vec d);
void VtoM (Vec a, Vec b,Vec c,Mat *M);
void Mprint (Mat *m);
void MmulM (Mat *p,Mat *q,Mat *R);
void Minv (Mat *m, Mat *W, float d);
float Mdet (Mat *m);
void MmulV (Mat *m, Vec d,Vec *e);
void VmulM (Mat *m, Vec d, Vec *e); 
int line2tri (Vec a, Vec b, Vec c, Vec d, Vec e); 
int line2line (Vec a, Vec b, Vec c, Vec d,float s);
float lineOline (Vec a, Vec b, Vec c,Vec d, Vec *box);
int line2dot (Vec a, Vec b,Vec c, float s); 
float dotOline (Vec a, Vec b, Vec c, Vec *e);
float torsion (Vec a, Vec b, Vec c, Vec d);
float	angle1pi (float s,float c);
float	angle2pi (float s, float c);
float angdif (float a, float b);

#endif

