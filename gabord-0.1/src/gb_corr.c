#include "mpp.h"
//#include <param.h>
/*
 * Calculate the inner product of two given indeces for two waveforms
 *
 * inputs:
 *	index1:	the index for the first waveform (INDEX)
 *	index2:	the index for the second waveform (INDEX)
 *	nL:	the maximum level for the transformation (int)
 *	nS:	the subsampling level (int)
 *	nDeltaP:the delta p (int)
 *	nDeltaW:the delta w (int)
 *	nBound:the bound for the index (int)
 *	pnAep	the array for the length of k for each n (int *)
 *	pnIep	the array for the index of n (int *)
 *	pfCE1:  the array that contains the cos and sin on integers (double *)
 *	pfCE2:  the array that contains the cos and sin on
 *		k2^n/(1+2^(2n)) (double *)
 *	pfCE3:	the array that contains the cos and sin on
 *		k/(1+2^(2n)) (double *)
 *	pfG:	the array that contains the Gaussian on
 *		k/(2^nS sqrt(1+2^(2n)))
 * output:
 *	pfReal:	the real part of the inner product
 *	pfImaginary:	the imaginary part of the inner product
 */
double GaborGetInnerProd(index1,index2,nL,nST,nSF,nDeltaP,nDeltaW,
	nBound,C,pnAep,
	pnIep,pfCE1,pfCE2,pfCE3,pfG,pfReal,pfImaginary)
INDEX index1;
INDEX index2;
int nL;
int nST, nSF;
int nDeltaP, nDeltaW;
int nBound;
double C;
int *pnAep;
int *pnIep;
double *pfCE1,*pfCE2,*pfCE3,*pfG;
double *pfReal, *pfImaginary;
{
    int m, q;
    int j1, j2, p1, p2;
    int kt, nt, kw, nw, kw1, wk, wk1, tk, twk;
    int l1, l2, h1, h2;
    int abskt, abskw;
    int alpha, beta, gamma;
    int  N;
    int nStepM, nStepQ;
    int nIndex;
    int maxlh;
    double tmpt, tmpw, tmptw, tmpi, tmpr;
    double fCos, fSin;
    void GaborGetIndexForGauss(), GaborGetIndexForCExp();
    double GaborGetGauss();
    void GaborGetCExp();

    tmpr=tmpi=0.0;
    j1 = (int)index1->octave;
    j2 = (int)index2->octave;
    p1 = (int)index1->position;
    p2 = (int)index2->position;
    N = 1<<nL;
    maxlh = MAX(nST,nSF);
    l1 = j1>nST ? nST:j1;
    l2 = j2>nST ? nST:j2;
    h1 = nL-j1>nSF ? nSF:nL-j1;
    h2 = nL-j2>nSF ? nSF:nL-j2;
    if (j2>j1)
	nStepM = 1<<(nL-j1+maxlh);
    else
	nStepM = 1<<(nL-j2+maxlh);
    if (j2<j1)
	nStepQ = 1<<(j1+maxlh);
    else
	nStepQ = 1<<(j2+maxlh);

    GaborGetIndexForGauss(j1,j2,nDeltaP,nL,maxlh,l1,l2,
				-nBound,&kt,&nt);
    GaborGetIndexForGauss(nL-j1,nL-j2,nDeltaW,nL,maxlh,h1,h2,
				-nBound,&kw1,&nw);
    wk1 = nDeltaW - nBound*N;
    tk = nDeltaP - nBound*N;
    for (m=(-nBound);m<=nBound;m++)
	{
	abskt = abs(kt);
	if (abskt>=pnAep[nt])
	    {
	    kt += nStepM;
	    tk += N;
	    continue;
	    }
	tmpt = GaborGetGauss(pfG,pnIep,abskt,nt);
	kw = kw1;
	wk = wk1;
	for (q=(-nBound);q<=nBound;q++)
	    {
	    abskw = abs(kw);
	    if (abskw>=pnAep[nw])
		{
		kw += nStepQ;
		wk += N;
		continue;
		}
	    tmpw = GaborGetGauss(pfG,pnIep,abskw,nw);
	    tmptw = tmpw*tmpt;
	    twk = wk*tk;
	    GaborGetIndexForCExp(nt,abs(twk),&alpha,&beta,&gamma);
	    GaborGetCExp(alpha,beta,gamma,N,nt,
		pfCE1,pfCE2,pfCE3,&fCos,&fSin);
	    tmpr += fCos*tmptw;
	    if (twk<0)
	        tmpi -= fSin*tmptw;
	    else
	        tmpi += fSin*tmptw;
	    kw += nStepQ;
	    wk += N;
	    }
	kt += nStepM;
	tk += N;
	}
    if (j2<j1)
	{
	nIndex = p2*abs(nDeltaW);
	nIndex %= N;
	fCos = pfCE1[nIndex];
	if (nDeltaW>0)
	    fSin = pfCE1[nIndex+N];
	else
	    fSin = -pfCE1[nIndex+N];
	*pfReal = C*(tmpr*fCos+tmpi*fSin);
	*pfImaginary = C*(-tmpr*fSin+tmpi*fCos);
	}
    else
	{
	nIndex = p1*abs(nDeltaW);
	nIndex %= N;
	fCos = pfCE1[nIndex];
	if (nDeltaW>0)
	    fSin = pfCE1[nIndex+N];
	else
	    fSin = -pfCE1[nIndex+N];
	*pfReal = C*(tmpr*fCos-tmpi*fSin);
	*pfImaginary = -C*(tmpr*fSin+tmpi*fCos);
	}
}
/*
 * get the index for the look up table for the inner product
 *
 * inputs:
 *      j1; the octave for the first waveform (int)
 *      j2; the octave for the second waveform (int)
 *      nDeltaP: the difference of the two positions of the two
 *		 waveforms (int)
 *      m:  the multiple of N (int)
 *      nL: the maximum level for the transformation (int)
 *      nS: the level that starts the subsampling
 *
 * output:
 *      (k,n)   the index for the lookup table (int, int)
 */
void GaborGetIndexForGauss(j1,j2,nDeltaP,nL,maxlh,l1,l2,m,k,n)
int j1, j2;
int nDeltaP;
int nL;
int maxlh, l1, l2;
int m;
int *k, *n;
{
    int tmp;
 
    *n = abs(j2-j1);
    tmp = abs(nDeltaP);
    if (j2>j1)
        {
        tmp >>= j1-l1;
        if (nDeltaP>0)
            *k =  tmp+m*(1<<(nL-j1+l1));
        else
            *k =  -tmp+m*(1<<(nL-j1+l1));
        *k <<= maxlh-l1;
        }
    else
        {
        tmp >>= j2-l2;
        if (nDeltaP>0)
            *k =  tmp+m*(1<<(nL-j2+l2));
        else
            *k =  -tmp+m*(1<<(nL-j2+l2));
        *k <<= maxlh-l2;
        }
}
/*
 * Get the Gaussian value
 *
 * input:
 *	t:	the t square (double) (will be remaved)
 *	(k,n)   the index for the look up table (int, int)
 *
 * output:
 *	g:	the corresponding Gaussian g(sqrt(t))
 */
double GaborGetGauss(pfG,pnIep,k,n)
double *pfG;
int *pnIep;
int k;
int n;
{
     return(pfG[pnIep[n]+k]);
}
/*
 * get the array for the length of k for each n
 *
 * Inputs:
 *	nL:	the maximum number of levels in transformation
 *		signal size is 2^nL
 *	maxlh:	max(h,l), where l is the subsampling level for time and
 *		h is the subsamling level for frequency
 *	epslon:	the given precision (double)
 *
 * Output:
 *	pnAep:	the corresponding array (int *)
 *	pnIep:	the index for the array (int *)
 *
 * Bugs:
 * 	pnAep and pnIep need to be allocated with size of nL-1
 */
void GaborGetNAep(nL,maxlh,epslon,pnAep,pnIep)
int nL;
int maxlh;
double epslon;
int *pnAep, *pnIep;
{
    int i, sum=0;

    if (pnAep == (int *)NULL || pnIep == (int *)NULL)
	perror("Null inputs!");
    pnIep[0] = 0;
    for (i=0;i<nL-1;i++)
	{
	pnAep[i] = (int)(sqrt((double)(1+(1<<(2*i)))*(-log(epslon))/M_PI)
		+1.0)*(1<<maxlh);
	sum += pnAep[i];
	pnIep[i+1] = sum;
	}
}
/*
 * get the index for the complex exponential
 */
void GaborGetIndexForCExp(nDeltaJ,tw,alpha,beta,gamma)
int nDeltaJ;
int tw;
int *alpha, *beta, *gamma;
{
    int tmp1, tmp2, k;

    k = tw;
    tmp2 = 1+(1<<(2*nDeltaJ));
    *alpha = k/tmp2;
    tmp1 = k%tmp2;
    tmp2 = 1<<nDeltaJ;
    *beta = tmp1/tmp2;
    *gamma = tmp1%tmp2;
}
/*
 * get the cos
 */
void GaborGetCExp(alpha,beta,gamma,N,n,pfCE1,pfCE2,pfCE3,pfCos,pfSin)
int alpha, beta, gamma;
int N, n;
double *pfCE1, *pfCE2, *pfCE3;
double *pfCos, *pfSin;
{
    int index1, index2;
    int size;
    double a1, a2, a3, b1, b2, b3;
    double tmp1, tmp2;
    //double value;

    index1 = (1<<(n+1))-2;
    index2 = index1 + (n<<1);
    size = 1<<n;
    alpha %= N;
    a1 = pfCE1[alpha];
    b1 = pfCE1[alpha+N];
    a2 = pfCE2[index2+beta];
    b2 = pfCE2[index2+size+1+beta];
    a3 = pfCE3[index1+gamma];
    b3 = pfCE3[index1+size+gamma];

    tmp1 = a1*a2-b1*b2;
    tmp2 = a1*b2+a2*b1;
    *pfCos = a3*tmp1-b3*tmp2;
    *pfSin = a3*tmp2+b3*tmp1;
}
/*
 * get the sin
 */
double GaborGetSin(alpha,beta,gamma,N,n,pfCE1,pfCE2,pfCE3)
int alpha, beta, gamma;
int N, n;
double *pfCE1, *pfCE2, *pfCE3;
{
    int index1, index2;
    double value;

    index1 = (1<<(n+1))-2;
    index2 = index1 + (n<<1);
    value = pfCE1[alpha%N+N]*pfCE2[index2+beta]*pfCE3[index1+gamma];

    return(value);
}
/*
 * get the gaussina array
 *
 * input:
 *	nL: 2^nL is the signal size
 *	npAep:	the array for the k bound for each n (int *)
 *	length: the length of the array (int)
 *	maxlh:	max(l,h), where l is the subsample level for time and
 *		h is the subsample level for frequency
 *
 * output:
 *	pfG: the array of the gaussian values at k/sqrt(1+2^(2n)) (double *)
 */
double *GaborGetGaussianArray(nL,pnAep,length,maxlh)
int nL;
int *pnAep;
int length;
int maxlh;
{
    double *pfG;
    int k, n;
    double d, t, *v;

    pfG = (double *)malloc(sizeof(double)*length);
    v = pfG;
    for (n=0;n<nL-1;n++)
	{
	d = (double)((1+(1<<(n<<1)))*(1<<(maxlh<<1)));
	for (k=0;k<pnAep[n];k++)
	    {
	    t = (double)k*(double)k/d;
	    *v++ = exp(-M_PI*t);
	    }
	}

    return(pfG);
}
/*
 * get the complex exponential array (the first lattice)
 *
 * input:
 *	nL: 2^nL is the size of the signal (int)
 *
 * output:
 *	pfCE1: the complex exponential values at k 2^n/(1+2^(2*n)) (double *)
 */
double *GaborGetCE1(nL)
int nL;
{
    double *pfCE1;
    int n, N;
    double t, *vi, *vr, c;

    N = 1<<nL;
    pfCE1 = (double *)malloc(sizeof(double)*(N<<1));
    if (pfCE1 == (double *)NULL)
	perror("mem. alloc. failed!");
    c = 2.0*M_PI/N;
    vr = pfCE1;
    vi = pfCE1+N;
    for (n=0;n<N;n++)
	{
	t = (double)n;
	*vr++ = cos((double)(c*t));
	*vi++ = sin((double)(c*t));
	}

    return(pfCE1);
}
/*
 * get the complex exponential array (the second lattice)
 *
 * input:
 *	nL: 2^nL is the size of the signal (int)
 *
 * output:
 *	pfCE2: the complex exponential values at k 2^n/(1+2^(2*n)) (double *)
 */
double *GaborGetCE2(nL)
int nL;
{
    double *pfCE2;
    int length, d1, n, k, N;
    double t, d, *vi, *vr, c;

    length = (1<<nL)-2+((nL-1)<<1);
    pfCE2 = (double *)malloc(sizeof(double)*length);
    if (pfCE2 == (double *)NULL)
	perror("mem. alloc. failed!");
    N = 1<<nL;
    c = 2.0*M_PI/N;
    vr = pfCE2;
    for (n=0;n<nL-1;n++)
	{
	d1 = 1<<n;
	d = (double)(1+(d1<<n));
	vi = vr+d1+1;
	for (k=0;k<=d1;k++)
	    {
	    t = (double)(k*d1)/d;
	    *vr++ = cos((double)(c*t));
	    *vi++ = sin((double)(c*t));
	    }
	vr = vi;
	}

    return(pfCE2);
}
/*
 * get the complex exponential array (the third lattice)
 *
 * input:
 *	nL: 2^nL is the size of the signal (int)
 *
 * output:
 *	pfCE3: the complex exponential values at k 2^n/(1+2^(2*n)) (double *)
 */
double *GaborGetCE3(nL)
int nL;
{
    double *pfCE3;
    int length, d1, n, k, N;
    double t, d, *vi, *vr, c;

    length = (1<<nL)-2;
    pfCE3 = (double *)malloc(sizeof(double)*length);
    if (pfCE3 == (double *)NULL)
	perror("mem. alloc. failed!");
    N = 1<<nL;
    c = 2.0*M_PI/N;
    vr = pfCE3;
    for (n=0;n<nL-1;n++)
	{
	d1 = 1<<n;
	d = (double)(1+(d1<<n));
	vi = vr+d1;
	for (k=0;k<d1;k++)
	    {
	    t = (double)k/d;
	    *vr++ = cos((double)(c*t));
	    *vi++ = sin((double)(c*t));
	    }
	vr = vi;
	}

    return(pfCE3);
}
/*
 * get the coeff array (sqrt(2^n/((1+2^(2n))2 M_PI)
 *
 * input:
 *	nL: 2^nL is the signal size (int)
 *
 * output:
 *	pfC:	the corresponding array (double *)
 */
double *GaborGetCArray(nL)
int nL;
{
    double *pfC;
    int n;

    pfC = (double *)malloc(sizeof(double)*(nL-1));
    if (pfC == (double *)NULL)
	perror("mem. alloc. failed!");

    for (n=0;n<nL-1;n++)
	pfC[n] = (double)sqrt((double)(1<<n)/(double)((1+(1<<(n<<1)))*
		2.0*M_PI));

    return(pfC);
}
/*
 * get the constants for the inner products
 *
 * inputs:
 *	pfGNorm:	the array containing the discrete norm of the
 *		gaussian (double *)
 *	j1:	the octave for the first waveform (int)
 *	j2:	the octave for the second waveform (int)
 *	pfC:	the array containing sqrt(2^n/((2^(2n)+1)2 M_PI)) (double *)
 *
 * output:
 *	C:	the corresponding coeff. (double)
 */
double GaborGetCoeff(pfGNorm,j1,j2,pfC)
double *pfGNorm;
int j1, j2;
double *pfC;
{
    double k1, k2, C;
    int n;
    //double A;

    k1 = pfGNorm[j1];
    k2 = pfGNorm[j2];
    n = abs(j2-j1);
/*
    A = (double)((1<<(j1<<1))+(1<<(j2<<1)));

    C = k1*k2*(double)sqrt((double)(1<<(j1+j2))/(double)(A*2.0*M_PI));
*/
    C = k1*k2*pfC[n];

    return(C);
}
/*
 * get the index bound for the sum
 *
 * inputs:
 *	epslon:	the precision (double)
 *
 * output:
 *	nBound:	the bound (int)
 */
int GaborGetBound(epslon)
double epslon;
{
    int nBound;

    nBound = (int)sqrt(-log((double)epslon)/M_PI)+1;

    return(nBound);
}
/*
 * get the bound coefficient
 */
double *GaborGetB(nL,epslon)
int nL;
double epslon;
{
    double *bound, c;
    int i;

    bound = (double *)malloc(sizeof(double)*(nL-1));
    if (bound == (double *)NULL)
	perror("mem. alloc. failed!");

    c = -log(epslon)/M_PI;
    for (i=0;i<nL-1;i++)
	bound[i] = (double)sqrt((double)(1+(1<<(i<<1)))*c);

    return(bound);
}
/*
 * end of ng_corr.c
 */
