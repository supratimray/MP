/*********************************************************
*
* SRFFT.C  - Split-Radix Fast Fourier Transform
* 
* This is a decimation in frequency version as described
* by Malvar in his book on lapped orthogonal transforms.
* 
* Author: Henrique Malvar
* Date: 10/8/1991
*
* Usage:  srfft(xr, xi, logm) -- for direct DFT
*         srifft(xr,xi,logm) -- for inverse DFT
*
* Arguments:  xr(double)  input and output vector, length M,
*                real part  
*             xi(double) imaginary part
*              logm (int)  log(base 2) of vector length M
*
*
************************************************************/


#include <stdio.h>
#include <math.h>
#ifdef sun4
#include <alloca.h>
#endif
#include <stdlib.h>
#include "mpp.h"


#define MAXLOGM  25
#define TWOPI  6.28318530717958647692
#define SQHALF   0.707106781186547524401

static void error_exit();
void BR_permute();
void srrec();
void srfft();
void srifft();

/*
fft(args)
char *args;
{
SIGNAL sig_in_com, sig_out_com, sig_in, sig_out, out_com, out;
int size, logm, typ, newsize, typ2, newsize2, i, tvar;
args= check_format(args, "%I %I %O %O %?d1 %?d0", &sig_in, &sig_in_com, &out, &out_com, &typ, &typ2);

tvar = 0;

while (args = strchr(args, '-'))
    switch (*++args) {
    case 'z':
      tvar = 1;
      fprintf( foutput, "FFT begins at point 0 \n");
      break;
    default:
      error_option(*args);
    }

size = sig_in->size;
if (sig_in_com->size != size)
  error("You idiot! Real and complex signals must have the same sizes. What do you expect me to do?");
logm = find2power(sig_in->size);

newsize = ( 1 << logm);
sig_out_com = new_struct_signal();
sig_out = new_struct_signal();
change_signal(sig_out_com, newsize);
change_signal(sig_out, newsize);

for (i=0;i<size;i++){
if (tvar == 0){
  if (i >= size/2){
    sig_out->values[i - size/2] =  sig_in->values[i];
    sig_out_com->values[i - size/2] = sig_in_com->values[i];
  }
  else {
	   sig_out->values[i + newsize - size/2] =  sig_in->values[i];
	   sig_out_com->values[i + newsize - size/2] = sig_in_com->values[i];
	     } 
}
  else if (tvar == 1){
    sig_out->values[i] =  sig_in->values[i];
    sig_out_com->values[i] = sig_in_com->values[i];
  }
}

if(newsize > size){
    switch(typ){
      case 0:
	for(i=size/2; i < newsize - size/2 ; i++){
	  sig_out->values[i] = 0.0; 
	  sig_out_com->values[i] = 0.0; 
	}
        break;
      case 1:
        for(i=size/2; i < newsize - size/2; i++){
           sig_out->values[i] = sig_out->values[i-1] - (i-size/2+1)*(sig_out->values[newsize - size/2]-sig_out->values[size/2 -1])/(newsize - size + 1);
           sig_out_com->values[i] = sig_out_com->values[i-1] - (i-size/2+1)*(sig_out_com->values[newsize - size/2] - sig_out_com->values[size/2 - 1])/(newsize - size + 1);
        	}			
        break;
  default:
        newsize2 = newsize - size;
        newsize2 /= 2;
        for(i=size/2; i < size/2 + newsize2; i++){
  	     sig_out->values[i] = sig_out->values[size - i - 1];
	     sig_out_com->values[i] = sig_out_com->values[size - i - 1];
         	}	
        for(i=size/2 + newsize2; i < newsize - size/2; i++){
  	     sig_out->values[i] = sig_out->values[2*newsize - size + i];
	     sig_out_com->values[i] = sig_out_com->values[2*newsize -size + i];
         	}	
         break;
      }
  }

if (typ2 == 0)
srfft(sig_out->values, sig_out_com->values, logm);
else if (typ2 == 1)
srifft(sig_out->values, sig_out_com->values, logm);

change_signal(out_com, newsize);
change_signal(out, newsize);

for (i=0;i<newsize;i++){
if(tvar == 0){
  if (i >= newsize/2){
    out->values[i - newsize/2] =  sig_out->values[i];
    out_com->values[i - newsize/2] = sig_out_com->values[i];
  }
  else {
out->values[i + newsize/2] =  sig_out->values[i];
out_com->values[i + newsize/2] = sig_out_com->values[i];
	     }
}
else if(tvar == 1){
    out->values[i] =  sig_out->values[i];
    out_com->values[i] = sig_out_com->values[i];
  }
}
delete_signal(sig_out_com);
delete_signal(sig_out);
return;
}
*/
/*********************************************************
* FFT(input signal, &output signal)                   
**********************************************************/


FFT(sigin, psigout)
SIGNAL sigin, *psigout;
{
  double *xr, *xi;
  int logm,  HalfSize;

  if(sigin == NULL) error("FFT: sigin is NULL \n");
  if(sigin->size %2 != 0) 
    error("FFT: size of complex signal should be even \n");
  HalfSize = sigin->size/2;
  logm = find2power(sigin->size);
  if(*psigout == NULL) 
    *psigout = new_struct_signal();
  if(sigin != *psigout) 
    sig_copy(sigin, *psigout);
  xr = (*psigout)->values;
  xi = (*psigout)->values + HalfSize;
  srfft(xr,xi,logm-1);
/*  (*psigout)->values = xr; */
/*  (*psigout)->values + HalfSize = xi; */
}


/*********************************************************
* IFFT(input signal, &output signal)                   
**********************************************************/

IFFT(sigin, psigout)
     SIGNAL sigin, *psigout;
{
  double *xr, *xi;
  int logm, find2power(), HalfSize;

  if(sigin == NULL) error("FFT: sigin is NULL \n");
  if(sigin->size %2 != 0) 
    error("FFT: size of complex signal should be even \n");
  HalfSize = sigin->size/2;
  logm = find2power(sigin->size);
  if(*psigout == NULL) 
    *psigout = new_struct_signal();
  
  if(sigin != *psigout)    
    sig_copy(sigin, *psigout);

  xr = (*psigout)->values;
  xi = (*psigout)->values + HalfSize;
  srifft(xr,xi,logm-1);
}
  
  
  

/**********************************************************
* Error exit for program abortion
**********************************************************/

static void error_exit()
{
  exit(1);
}

/**********************************************************
* Data unshuffling according to bit-reversed indexing.
*
* Bit reversal is done using Evans' algorithm. (Ref: DMW Evans,
* "An Improved Digit-Reversal Permutation Algorithm",
*  IEEE Trans. ASSP, Aug 1987, pp. 1120-25.
*
***********************************************************/

static int brseed[256]; /* Evans' seed table */
static int brsflg;   /* flag for table building */

void BR_permute(x, logm)
double *x;
int logm;
{
register int lg2, n;
int i,j,imax;
int off,fj, gno, *brp;
double tmp, *xp, *xq;

lg2 = logm >> 1;
n = 1 << lg2;

if (logm & 1) lg2++;

/* create seed table if not yet built */
if (brsflg != logm) {
  brsflg = logm;
  brseed[0] = 0;
  brseed[1] = 1;
  for(j=2; j <= lg2; j++){
    imax =  1 << (j - 1);
    for (i=0; i < imax; i++){
      brseed[i] <<= 1;
      brseed[i + imax] = brseed[i] + 1;
   }
 }
}

/* unshuffling loop */
for(off = 1; off < n; off++){
fj = n*brseed[off]; i = off; j = fj;
tmp = x[i]; x[i] = x[j]; x[j] = tmp;
xp = &x[i];
brp = &brseed[1];
for(gno = 1; gno < brseed[off]; gno++){
  xp += n;
  j = fj + *brp++;
  xq = x + j ;
  tmp = *xp; *xp = *xq; *xq = tmp;
}
}


}

/**************************************************************
*                                                             *
* Recursive part of split radix FFT algorithm
*
***************************************************************/

void srrec(xr, xi, logm)
double *xr, *xi;
int logm;
{

static int m, m2, m4, m8, nel, n;
static double *xr1, *xr2, *xi1, *xi2;
static double *cn, *spcn, *smcn, *c3n, *spc3n, *smc3n;
static double tmp1, tmp2, ang, c, s;
static double *tab[MAXLOGM];

/* check range of logm */

if ((logm < 0) || (logm > MAXLOGM)){
 fprintf( foutput, "Error: SRFFT logm = %d is out of bounds [%d %d]\n", logm, 0, MAXLOGM);
  error_exit();
}

/* compute trivial cases */
if (logm < 3){
if (logm == 2){
xr2 = xr + 2;
xi2= xi + 2;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
xr1 = xr + 1;
xi1 = xi + 1;
xr2++;
xi2++;
tmp1 = *xr1 + *xr2;
*xr2 = *xr1 - *xr2;
*xr1 = tmp1;
tmp1 = *xi1 +  *xi2;
*xi2 = *xi1 - *xi2;
*xi1 = tmp1;
xr2 = xr + 1;
xi2 = xi + 1;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
xr1 = xr +2;
xi1 = xi+2;
xr2 = xr+3;
xi2 = xi + 3;
tmp1 = *xr1 + *xi2;
tmp2 = *xi1 + *xr2;
*xi1 = *xi1 - *xr2;
*xr2 = *xr1 - *xi2;
*xr1 = tmp1;
*xi2 = tmp2;
return;
}
else if (logm == 1){
xr2 = xr + 1;
xi2 = xi + 1;
tmp1 = *xr + *xr2;
*xr2 = *xr - *xr2;
*xr = tmp1;
tmp1 = *xi + *xi2;
*xi2 = *xi - *xi2;
*xi = tmp1;
return;
}
else if (logm == 0) return;
}

/* compute a few constants */

m = 1 << logm; m2 = m / 2; m4 = m2 / 2; m8 = m4 / 2;


/* build tables of butterfly coefficients if necessary */
if ((logm >=4) && (tab[logm-4] == NULL)){
/* allocate memory for tables */
nel = m4 - 2;
if ((tab[logm-4] = (double *) calloc(6 * nel, sizeof(double)))
    == NULL){
  error_exit();
}



/* initialize pointers */
cn = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
c3n = smcn + nel; spc3n = c3n + nel; smc3n = spc3n + nel;


/* compute tables */
for (n = 1; n < m4; n++){
  if (n == m8) continue;
  ang = n * TWOPI / m;
  c = cos(ang); s = sin(ang);
  *cn++ = c; *spcn++ = -(s + c); *smcn++ = s - c;
  ang = 3*n*TWOPI/m;
  c = cos(ang); s = sin(ang);
  *c3n++ = c; *spc3n++ = -(s + c); *smc3n++ = s - c;
}
}


/* step 1 */
xr1 = xr; xr2 = xr1 + m2;
xi1 = xi; xi2 = xi1 + m2;

for (n =0; n < m2; n++){
  tmp1 = *xr1 + *xr2;
  *xr2 = *xr1 - *xr2;
  *xr1 = tmp1;
  tmp2 = *xi1 + *xi2;
  *xi2 = *xi1 - *xi2;
  *xi1 = tmp2;
  xr1++; xr2++; xi1++; xi2++;
} 

/* Step 2 */
xr1 = xr + m2;  xr2 = xr1 + m4;
xi1 = xi + m2;  xi2 = xi1+ m4;
for (n = 0; n < m4; n++){
  tmp1 = *xr1 + *xi2;
  tmp2 = *xi1 + *xr2;
  *xi1 = *xi1 - *xr2;
  *xr2 = *xr1 - *xi2;
  *xr1 = tmp1;
  *xi2 = tmp2;
  xr1++; xr2++; xi1++; xi2++;
}

/* Steps 3&4 */
xr1 = xr + m2; xr2 = xr1 + m4;
xi1 = xi + m2; xi2 = xi1 + m4;

if (logm >= 4) {
  nel = m4 -2;
  cn = tab[logm-4]; spcn = cn + nel; smcn = spcn + nel;
  c3n = smcn + nel; spc3n = c3n + nel; smc3n = spc3n + nel;
}

xr1++; xr2++; xi1++; xi2++;

for(n=1; n < m4; n++){
   if (n == m8){
     tmp1 = SQHALF*(*xr1 + *xi1);
     *xi1 = SQHALF*(*xi1 - *xr1);
     *xr1 = tmp1;
     tmp2 = SQHALF*(*xi2 - *xr2);
     *xi2 = -SQHALF*(*xr2 + *xi2);
     *xr2 = tmp2;
   }
   else {
     tmp2 = *cn++ *(*xr1 + *xi1);
     tmp1 = *spcn++ * *xr1 + tmp2;
     *xr1 = *smcn++ * *xi1 + tmp2;
     *xi1 = tmp1;
     tmp2 = *c3n++ * (*xr2 + *xi2);
     tmp1 = *spc3n++ * *xr2 + tmp2;
     *xr2 = *smc3n++ * *xi2 + tmp2;
     *xi2 = tmp1;
   }
xr1++; xr2++; xi1++; xi2++;
 }


/* call ssrec again with half DFT length */
srrec(xr, xi, logm -1);
/* call ssrec again twice with one quarter DFT length.
Constants have to be recomputed because they are static! */

m = 1 << logm; m2 = m/2;
srrec(xr+m2, xi+m2, logm - 2);

m = 1 << logm; m4 = 3*(m/4);
srrec(xr+m4, xi+m4, logm - 2);
}

/*****************************************************************
*  Direct transform                                               *
*
*******************************************************************/
	
void srfft(xr, xi, logm)
double *xr, *xi;
int logm;
{
/* call recursive routine */
srrec(xr, xi, logm);

/* output array unshuffling using bit-reversed indices */
if (logm > 1){
  BR_permute(xr, logm);
  BR_permute(xi, logm);
}
}

/**************************************************************
*
*Inverse transform. Uses Duhamel's trick (Ref: P. Duhamel
* et al. "On Computing the Inverse DFT", IEEE Trans. ASSP
* Feb. 1988, pp. 285-286.
*
****************************************************************/
void srifft(xr, xi, logm)
double *xr, *xi;
int logm;
{
  int i,m;
  double fac, *xrp, *xip;

/* call direct FFT, swapping real and imaginery addresses */
  srfft(xi, xr, logm);
  
/* Normalization */
  m = 1 << logm;
  fac = 1.0/m;
  xrp = xr; xip = xi;

for (i=0; i < m; i++){
  *xrp++ *= fac;
  *xip++ *= fac;
}
}
