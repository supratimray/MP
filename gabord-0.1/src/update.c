#include "mpp.h"
#define mqMax 4
#define mqMax2 2*mqMax
#define mqMax3 mqMax2+mqMax
#define mqMax4 4*mqMax
#define SSFACTOR(j,l)   (((j)>(l))?(l):(j))     /* subsample factor */
#define STEP(j,l)       (1<<((j)-(l)))  /* calcalate the step for $\delta$ */
#define AB(j1,j2)       ((1<<((j1)<<1))+(1<<((j2)<<1))) /* compute the A or B */
#define DELTAMAXP(j,N,a)      (MIN((N),((long)(a)<<(j))))
#define DELTAMAXN(j,N,a)      (MAX(-(N),-((long)(a)<<(j))))
/* real part of complex multiplication */
#define CMREAL(a1,b1,a2,b2)	((a1)*(a2)-(b1)*(b2))
/* imaginary part of complex multiplication */
#define CMIMAGINARY(a1,b1,a2,b2)	((a1)*(b2)+(a2)*(b1))
/* positive index for the real part */
#define POSINDREAL		mqMax+1
/* positive index for the imaginary part */
#define POSINDIMAG		mqMax3+1
/* negative index for the real part */
#define NEGINDREAL		mqMax-1
/* negative index for the imaginary part */
#define NEGINDIMAG		mqMax3-1
/* positive zero index for the real part */
#define POSZEROREAL		mqMax
/* positive zero  index for the imaginary part */
#define POSZEROIMAG		mqMax3
/* negative index for the real part */
#define NEGZEROREAL		mqMax
/* negative index for the imaginary part */
#define NEGZEROIMAG		mqMax3
/* the left bound of m or q */
#define LB(L,j,a,d)     (-((((long)(a)<<(j))+(d))>>(L)))
#define RB(L,j,a,d)   (((((long)(a)<<(j))-(d))>>(L)))
/* calculate the index for the Gaussian */
#define GINDX(a,b,c,d)   ((((a)<<(b))>>(c))+(d))
/* calculate the index for the transform given by $k_2$ and $p_2$ */
#define TINDXREAL(k2,p2,Ljl2)        (((k2)<<((Ljl2)+1))+(p2))
#define TINDXIMAG(k2,p2,Ljl2)         (((k2)<<((Ljl2)+1))+(1<<(Ljl2))+(p2))
/* max p2 and k2  */
#define P2MAX(L,j2,l2)  (1<<((L)-(j2)+(l2)))
#define K2MAX(j2,h2)    (1<<((j2)+(h2)-1))
/*
 * routine for update the Gabor transformation
 *
 * Inputs:
 *      trans	gabor transform (SIGNAL *)
 *      wrd    the wrd selected (WRD)
 *      L       $2^L$ is the signal size (int)
 *      l       subsampling octave for time domain (int)
 *      h       subsampling octave for the frequency domain (int)
 *      mqMax   the maximum number for m and q (global constant)
 *      mqMax2	2*mqMax
 *      $pfG[j_1][j_2][k] = e^{-\pi\frac{k}{2^(2j_1)+2^(2J_2)}$ (global variable)
 *      $pfC[n] = \sqrt(\frac{2^n}{1+2^{2n}}) $ (global variable)
 *      pfGNorm[j]	the discrete normalization factor of $g(t/2^j)$ (global variable)
 *      pfB[n]  $\sqrt (\frac{(1+2^{2n})log1/\epsilon}{\pi})$
 *	$pfCE1[n] = e^{-i\frac{2\pi}{N}n}$
 *	pnIep:	the address of n for pfG (long *)
 *
 * Output:
 *      trans   updated gabor transform (SIGNAL *)
 *
 * Bugs:
 *	Assuming that the first wave form (j1,k1,p1) corresponding to
 * 	the fine grid and the second wave form (j2,k2,p2) corresponding to
 *	the coarse grid
 */
void UpdateGabor(trans,wrd, MinOctave, MaxOctave, L,l,h,
		lf, hf,pfG,pfC,pfFilterNorm,
		pfFilter,pfB,pfCE1,pnIep)
SIGNAL *trans;
WRD wrd;
int MinOctave, MaxOctave; /* Min and Max octave for the decomposition on Gabor functions */
int L;
int l, h;	/* oversampling octaves on the coarse grid */
int lf, hf;	/* oversampling octaves on the fine grid */
double *pfG;
double *pfC;
double *pfFilterNorm;
FILTER *pfFilter;
double *pfB;
double *pfCE1;
int *pnIep; //back to int for 64bits
{
/*
 * local variables for the complex exponential
 */
    static double f1_0[2];   /* $f_1(\delta,0,0)$ where $\delta=\delta k_0$
                              * or $\delta p_0$
                              */
    static double f2_0[2];   /* $f_2(\delta,\delta,0)$ where $\delta=\delta p_0$
			      * or $\delta k_0$
                              */
    static double fs2_0[2];  /* $f_2(step,\delta,0)$ */
    static double f3_0[2];   /* $f_3(step,0)$ where $step=stepP$ or $stepK$ */
/*
 * local variables for parameters
 */
    long N=1<<L;        /* signal size */
    int j1, j2;         /* octaves */
    int Lj1, Lj2;         /* $L-j_1$ and $L-j_2$ */
    long k2, p2;
    long k1, p1;/* the indeces for frequncy and time */
    long A, B;          /* the A, B */
    long dk0, dp0;      /* $\delta k_0$ and $\delta p_0$ */
    long dkMaxP, dpMaxP;/* bounds for $\delta k\geq 0$ and $\delta p\geq 0$ */
    long dkMaxN, dpMaxN;/* bounds for $\delta k\leq 0$ and $\delta p\leq 0$ */
    int h1, h2, l1, l2;  /* subsampling factors for frequency and time */
    long stepK, stepP;	/* the steps for dk and dp */
    double coeff;       /* coefficient */
    double phi;         /* $\phi$ */
    double cos_phi;     /* $cos(\phi)$ */
    double sin_phi;     /* $sin(\phi)$ */
    double cos_2phi;	/* $cos(2\phi)$ */
    double sin_2phi;	/* $sin(2\phi)$ */
    double sqrtN;
    int n;
    short dkFirst;
/*
 * temperary variables
 */
    INDEX index;
    /* $\delta k=0$ or $stepK/2$, $\delta p=0$ or $stepP/2$ */
    short dkZero, dpZero;
    long ntmp, k2Max, p2Max, kIndxStep, kIndx, pIndxStep, pIndx;
    double dtmp1, dtmp2, dReal1, dReal2, dImag1, dImag2, dCos, dSin, dCos2, dSin2;
    double *vr, *vi, *vg;
    double C; /* the normalization factors */
    double cpi; /* $2\pi 2^{2j_1}/A$ */
    unsigned long flag=0L;
    int minj, minLj, maxlh;
/*
 * functions used
 */
    double GaborGetCoeff();
    void genericLoop();
    void initValues();

    sqrtN = sqrt((double)N);
/*
 * get the information for the selected wrd
 */
    index = wrd->index;
    coeff = wrd->coeff*wrd->value;
    phi = index->phase;
    p1 = (long)index->position;
    k1 = (long)index->id;
    j1 = (int)index->octave;
    Lj1 = L-j1;
    l1 = SSFACTOR(j1,lf);
    h1 = SSFACTOR(Lj1,hf);
    maxlh = MAX(MAX(l,h),MAX(lf,hf));
/*
 * case 1: the selected wrd is dirac
 */
    if (j1==0)
        {
        if (fabs(phi)>1.0e-5) {
           coeff = -coeff;
        }
        for ( j2 = MinOctave; j2 <= MaxOctave; j2++) {
           Lj2 = L-j2;
           l2 = SSFACTOR(j2,l);
           h2 = SSFACTOR(Lj2,h);
           p2Max = P2MAX(L,j2,l2);
           k2Max = K2MAX(j2,h2);
           vg = pfFilter[j2]->values;
           vr = trans[j2 - MinOctave + 1]->values;
           kIndx = 0;
           kIndxStep = p1<<(Lj2-h2);
           for (k2=0; k2<=k2Max; k2++) {
              vi = vr+p2Max;
              ntmp = kIndx%N;
              cos_phi = coeff*pfCE1[ntmp];
              sin_phi = -coeff*pfCE1[N+ntmp];
              kIndx += kIndxStep;
              pIndx = p1+N;
              pIndxStep = 1<<(j2-l2);
              for (p2=0; p2<p2Max; p2++) {
                ntmp = pIndx%N;
		if (ntmp<pfFilter[j2]->size)
		    {
                    *vr++ -= vg[ntmp]*cos_phi;
                    *vi++ -= vg[ntmp]*sin_phi;
		    }
		else if (N-ntmp<pfFilter[j2]->size)
		    {
                    *vr++ -= vg[N-ntmp]*cos_phi;
                    *vi++ -= vg[N-ntmp]*sin_phi;
		    }
		else
		    {
                    vr++;
                    vi++;
		    }
                 pIndx -= pIndxStep;
              } /* end p2 loop */
              vr = vi;
           } /* end of k2 loop */
        } /* end of j2 loop */
        return;
        } /* end of j1==0 */
/*
 * case 2: the selected wrd is Fourier
 */
    cos_phi = coeff*cos(phi)/2.0;
    sin_phi = coeff*sin(phi)/2.0;
    if (j1==L)
	{
	for ( j2= MinOctave;j2 <= MaxOctave; j2++)
	    {
            Lj2 = L-j2;
            l2 = SSFACTOR(j2,l);
            h2 = SSFACTOR(Lj2,h);
	    p2Max = P2MAX(L,j2,l2);
	    k2Max = K2MAX(j2,h2);
	    vr = trans[j2 - MinOctave +1]->values;
	    vg = pfFilter[L-j2]->values;
            kIndxStep = Lj2-h2;
            pIndxStep = j2-l2;
	    for (kIndx=0;kIndx<=k2Max;kIndx++)
		{
		vi = vr+p2Max;
                k2 = kIndx<<kIndxStep;
		if (k1>=k2)
		    ntmp = k1-k2;
		else
		    ntmp = N-(k2-k1);
		if (ntmp<pfFilter[L-j2]->size)
		    dtmp1 = vg[ntmp]*sqrtN;
		else if (N-ntmp<pfFilter[L-j2]->size)
		    dtmp1 = vg[N-ntmp]*sqrtN;
		else
		    dtmp1 = 0.0;
/*
		if (k1>=k2)
		    dtmp1 = vg[k1-k2];
		else
		    dtmp1 = vg[N-(k2-k1)];
*/
		ntmp = N-k1-k2;
		if (ntmp<pfFilter[L-j2]->size)
		    dtmp2 = vg[ntmp]*sqrtN;
		else if (N-ntmp<pfFilter[L-j2]->size)
		    dtmp2 = vg[N-ntmp]*sqrtN;
		else
		    dtmp2 = 0.0;
/*
		if (ntmp==N)
		    dtmp2 = vg[0];
		else
		    dtmp2 = vg[ntmp];
*/
		for (pIndx=0;pIndx<p2Max;pIndx++)
		    {
		    p2 = pIndx<<pIndxStep;
		    p1 = abs(p2*(k1-k2));
		    dCos = pfCE1[p1%N];
		    dSin = pfCE1[N+p1%N];
		    p1 = abs(p2*(N-k1-k2));
		    dCos2 = pfCE1[p1%N];
		    dSin2 = pfCE1[N+p1%N];
		    dReal1 = dCos*dtmp1;
		    if (k2>k1)
		        dImag1 = -dSin*dtmp1;
		    else
			dImag1 = dSin*dtmp1;
		    dReal2 = dCos2*dtmp2;
		    if (k2>N-k1)
			dImag2 = -dSin2*dtmp2;
		    else
			dImag2 = dSin2*dtmp2;
		    *vr++ -= cos_phi*(dReal1+dReal2)+
				sin_phi*(dImag2-dImag1);
		    *vi++ -= cos_phi*(dImag1+dImag2)+
				sin_phi*(dReal1-dReal2);
		    } /* end of pIndx loop */
		vr = vi;
		} /* end of kIndx loop */
	    } /* end of j2 loop */
	return;
	} /* end of the Fourier case */
/*
 * case 3: the selected wrd is Gabor
 */
/*
 * DEBUG
 */
#ifdef DEBUG0
fprintf( foutput, "L = %d\t\tl = %d\th = %d\tl1 = %d\th1 = %d\n",L,l,h,l1,h1);
fprintf( foutput, "j1 = %d\tk1 = %d\tp1 = %d\tphi = %le\n",j1,k1,p1,phi);
fprintf( foutput, "coeff = %le\n",coeff);
#endif
    cos_2phi = cos(2.0*phi);
    sin_2phi = sin(2.0*phi);
    for ( j2=MinOctave; j2<= MaxOctave; j2++)
        {
        /* compute the parameters */
        n = abs(j2-j1);
        Lj2 = L-j2;
	minj = MIN(j1,j2);
	minLj = MIN(Lj1,Lj2);
        l2 = SSFACTOR(j2,l);
        h2 = SSFACTOR(Lj2,h);
        stepP = STEP(j2,l2);
        stepK = STEP(Lj2,h2);
        A = AB(j1,j2);
        B = AB(Lj1,Lj2);
        dp0 = delta0(stepP,p1);
        dk0 = delta0(stepK,k1);
        dkZero = (dk0==0) || (dk0==(stepK>>1));
        dpZero = (dp0==0) || (dp0==(stepP>>1));
        ntmp = N>>1;
        dpMaxP = DELTAMAXP(minj,ntmp,pfB[n]);
        dkMaxP = DELTAMAXP(minLj,ntmp,pfB[n]);
        dpMaxN = DELTAMAXN(minj,ntmp,pfB[n]);
        dkMaxN = DELTAMAXN(minLj,ntmp,pfB[n]);
        if (!dkZero) {
           while (dk0>dkMaxN) {
              dk0 -= stepK;
           } /* endwhile */
           if (dk0!=dkMaxN) {
              dk0 += stepK;
           }
        } else {
	   dkMaxP = MAX(dkMaxP,abs(dkMaxN));
	}/* endif */
        if (!dpZero) {
           while (dp0>dpMaxN) {
              dp0 -= stepP;
           } /* endwhile */
           if (dp0!=dpMaxN) {
              dp0 += stepP;
           }
        } else {
	   dpMaxP = MAX(dpMaxP,abs(dpMaxN));
	}/* endif */
        /* condition on looping $\delta k$ first*/
        if (dkZero && dpZero)
           {
           flag |= A_FLAG|B_FLAG; /* install the flags for dkZero and dpZero */
           dkFirst = j2<= j1+l2-h2;
           }
        else if (dkZero && !dpZero)
           {
           flag |= B_FLAG;
           flag ^= B_FLAG; /* delete the flag for dpZero */
           flag |= A_FLAG; /* install the flag for dkZero */
           dkFirst = j2<=j1+l2-h2+1;
           }
        else if (!dkZero && dpZero)
           {
           flag |= A_FLAG;
           flag ^= A_FLAG; /* delete the flag for dkZero */
           flag |= B_FLAG; /* install flag for dpZero */
           dkFirst = j2<=j1+l2-h2-1;
           }
        else {
           flag |= A_FLAG|B_FLAG;
           flag ^= A_FLAG|B_FLAG; /* delete the flags for dkZero and dpZero */
           dkFirst = j2<= j1+l2-h2;
        } /* endif */
        /* calculations of initial values */
        cpi = 2.0*M_PI*(double)(1<<(j1<<1))/(double)A;
        if (dkFirst) {
           /* condition on looping m first */
           if(j2<=L-j1)
	       flag |= M_FLAG; /* looping m first */
           else
               {
               /*
		* turn on the M_FLAG bit so that the next operation will
		* turn it off
		*/
               flag |= M_FLAG; 
	       flag ^= M_FLAG; /* looping q first */
               }
	   flag |= K_FLAG;
	   /*
	    * initial values for Procedure 1
	    */
	   /* compute f1_0 */
	   C = GaborGetCoeff(pfFilterNorm,j1,j2,pfC);
	   /*
	    * calculate the initial values
	    */
	   initValues(f1_0,f2_0,fs2_0,f3_0,dk0,dp0,p1,stepK,stepP,N,
			cpi,pfCE1,C*cos_phi,C*sin_phi,K_FLAG);
#ifdef DEBUG0
	fprintf( foutput, "L b4_1 = %d\n",L);        
#endif
	   genericLoop(trans, MinOctave, MaxOctave, f1_0,f2_0,fs2_0,f3_0,dk0,dp0,dkMaxP,
		dpMaxP,stepK,stepP,j1,j2,k1,p1,l1,l2,h1,h2,maxlh,
		A,B,L,flag,
		pfG,pnIep,pfCE1,pfB,cos_2phi,sin_2phi
#ifdef DEBUG0
                ,C*cos_phi,C*sin_phi
#endif
		);
        } else {
           /* condition on looping m first */
           if(j2<=L-j1)
              {
               flag |= M_FLAG;  /* turn on the M_FLAG bit so that the next operation will turn it off */
	       flag ^= M_FLAG; /* looping q first, this q is m in the formula */
               }
           else
	       flag |= M_FLAG; /* looping m first, this m is q in the formula */
           flag |= K_FLAG;
           flag ^= K_FLAG;
	   /*
	    * initial values for Procedure 1
	    */
	   /* compute f1_0 */
	   C = GaborGetCoeff(pfFilterNorm,j1,j2,pfC);
	   initValues(f1_0,f2_0,fs2_0,f3_0,dk0,dp0,p1,stepK,stepP,N,
			cpi,pfCE1,C*cos_phi,C*sin_phi,P_FLAG);
#ifdef DEBUG0
	fprintf( foutput, "L b4_2 = %d\n",L);        
#endif
		   genericLoop(trans, MinOctave, MaxOctave, f1_0,f2_0,fs2_0,f3_0,dp0,dk0,dpMaxP,
		dkMaxP,stepP,stepK,j1,j2,k1,p1,l1,l2,h1,h2,maxlh,
		A,B,L,flag,
		pfG,pnIep,pfCE1,pfB,cos_2phi,sin_2phi
#ifdef DEBUG0
		,C*cos_phi,C*sin_phi
#endif
		);
            } /* endif */
#ifdef DEBUG0
fprintf( foutput, "C = %le\n",C);
#endif
        } /* end of j2 loop */
} /* end of updateGabor() */
/*
 * calculate the $\delta p_0$ or $\delta k_0$
 * formula (1) in Part III
 *
 * Inputs:
 *	stepP:	the step size of the grid (int)
 *      p1:	the position of the gabor function
 *
 * Output:
 *      dp0:	the corresponding $\delta_0$
 */
int delta0(stepP,p1)
int stepP;
long p1;
{
    int delta;

    delta = (p1/stepP)*stepP - p1;
    if (delta<0)
	delta += stepP;

    return(delta);
}
/*
 * generic recursive formula
 *
 * Inputs:
 *	f[]:	array needed to be calculated (double *)
 *	posIndxR:	positive index of real part
 *	posIndxI:	positive index of imaginary part
 *	negIndxR:	negative index of real part
 *	negIndxI:	negative index of imaginary part
 *	fm:	the multiplier (double); fm[posIndxR] is the real part for positive index
 *		and fm[posIndxI] is the imaginary part; fm[negIndxR] is the real part
 *		for the negative index and fm[negIndxI] is the imaginary part for the
 *		negative index
 * Outputs:
 *	f[posIndxR], f[posIndxI], f[negIndxR],f[negIndxI]
 */
static void gnrRcsArray(f,posIndxR,posIndxI,negIndxR,negIndxI,fm)
double f[];
long posIndxR, posIndxI;
long negIndxR, negIndxI;
double fm[];
{
   double dtmp;

    dtmp = CMREAL(f[posIndxR],f[posIndxI],fm[posIndxR],fm[posIndxI]);
    f[posIndxI] = CMIMAGINARY(f[posIndxR],f[posIndxI],fm[posIndxR],fm[posIndxI]);
    f[posIndxR] = dtmp;
    dtmp = CMREAL(f[negIndxR],f[negIndxI],fm[negIndxR],fm[negIndxI]);
    f[negIndxI] = CMIMAGINARY(f[negIndxR],f[negIndxI],fm[negIndxR],fm[negIndxI]);
    f[negIndxR] = dtmp;
}
/*
 * generic recursive formula
 *
 * Inputs:
 *	f[]:	array needed to be calculated (double *)
 *	posIndxR:	positive index of real part
 *	posIndxI:	positive index of imaginary part
 *	negIndxR:	negative index of real part
 *	negIndxI:	negative index of imaginary part
 *	fmR, fmI:	the multiplier (double); fmR is the real part and fmI is the imaginary part.
 * Outputs:
 *	f[posIndxR], f[posIndxI], f[negIndxR],f[negIndxI]
 */
static void gnrRcs(f,posIndxR,posIndxI,negIndxR,negIndxI,fmR,fmI)
double f[];
long posIndxR, posIndxI;
long negIndxR, negIndxI;
double fmR, fmI;
{
   double dtmp;

    dtmp = CMREAL(f[posIndxR],f[posIndxI],fmR,fmI);
    f[posIndxI] = CMIMAGINARY(f[posIndxR],f[posIndxI],fmR,fmI);
    f[posIndxR] = dtmp;
    dtmp = CMREAL(f[negIndxR],f[negIndxI],fmR,fmI);
    f[negIndxI] = CMIMAGINARY(f[negIndxR],f[negIndxI],fmR,fmI);
    f[negIndxR] = dtmp;
}
/*
 * generic recursive formula
 *
 * Inputs:
 *	f[]:	array needed to be calculated (double *)
 *	posIndxR:	positive index of real part
 *	posIndxI:	positive index of imaginary part
 *	negIndxR:	negative index of real part
 *	negIndxI:	negative index of imaginary part
 *	fmR, fmI:	the multiplier (double);
 *			fmR is the real part and fmI is the imaginary part.
 * Outputs:
 *	f[posIndxR], f[posIndxI], f[negIndxR],f[negIndxI]
 */
void genericRecursion(f,posIndxR,posIndxI,negIndxR,negIndxI,fmR,fmI)
double f[];
long posIndxR, posIndxI;
long negIndxR, negIndxI;
double fmR, fmI;
{
    static long posPrevR, posPrevI, negPrevR, negPrevI; /* previous indeces */
    posPrevR = posIndxR-1;
    posPrevI = posIndxI-1;
    negPrevR = negIndxR+1;
    negPrevI = negIndxI+1;
    f[posIndxR] = CMREAL(f[posPrevR],f[posPrevI],fmR,fmI);
    f[posIndxI] = CMIMAGINARY(f[posPrevR],f[posPrevI],fmR,fmI);
    f[negIndxR] = CMREAL(f[negPrevR],f[negPrevI],fmR,-fmI);
    f[negIndxI] = CMIMAGINARY(f[negPrevR],f[negPrevI],fmR,-fmI);
}
/*
 * generic loop
 *
 * inputs:
 *	trans:	the gabor transformation (SIGNAL *)
 * 	MinOctave	Minimum of the octave decomposition on gabor functions
 *	MaxOctave	Maximum of the octave decomposition on gabor functions
 *	f1_0:	initial value for f1 (double *)
 *	f2_0:	initial value for f2 (double *)
 *	fs2_0:	initial value for f2 (double *)
 *	f3_0:	initial value for f3 (double *)
 *	dk0:	generic $\delta k_0$ (long)
 *	dp0:	generic $\delta p_0$ (long)
 *	dkMaxP:	max value of generic positive $\delta k$ (long)
 *	dkMaxN:	max value of generic negative $\delta k$ (long)
 *	dpMaxP:	max value of generic positive $\delta p$ (long)
 *	dpMaxN:	max value of generic negative $\delta p$ (long)
 *	stepK:	step size for generic $\delta k$ (long)
 *	stepP:	step size for generic $\delta p$ (long)
 *	j1,j2:	octaves (int)
 *	k1,k2:	frequencies (long)
 *	p1,p2:	positions (long)
 *	l1, l2:	subsampling rates for time (int)
 *	h1, h2:	subsampling rate for frequency (int)
 *	A, B:	the constant A and B (long)
 *	N:	signal size (long)
 *      pfG:    the array that contains the Gaussian on 
 *		$\frac{k}{2^j\sqrt(1+2^{2n))}$ (double *)
 *      pnIep:  the array for the index of n in pfG (int *)
 *      cR:     the real part of the coefficient (double)
 *      cI:     the imaginary part of the coefficient (double)
 *
 * Outputs:
 *	trans:	the updated transformation
 */
void genericLoop(trans, MinOctave, MaxOctave, f1_0,f2_0,fs2_0,f3_0,
		dk0,dp0,dkMaxP,dpMaxP,stepK,stepP,
		j1,j2,k1,p1,l1,l2,h1,h2,maxlh,
		A,B,L,flag,
		pfG,pnIep,pfCE1,pfB,cos_2phi,sin_2phi
#ifdef DEBUG0
		,cR,cI
#endif
		)
SIGNAL *trans;
int	MinOctave, MaxOctave; /* Min and Max octave for the decomposition on Gabor functions */
double f1_0[2];
double f2_0[2];
double fs2_0[2];
double f3_0[2];
long dk0, dp0;
long dkMaxP;
long dpMaxP;
long stepK, stepP;
int j1,j2;
long k1;
long p1;
int l1, l2, h1, h2;
int maxlh;
long A, B;
int L;
unsigned long flag;
double *pfG;
int *pnIep; //back to int for 64bits
double *pfCE1;
double *pfB;
double cos_2phi, sin_2phi;
#ifdef DEBUG0
double cR, cI;
#endif
{
/*
 * variables for the complex exponential
 */
    /* $f_4(\delta,0)$ where $\delta=\delta k_0$ or $\delta p_0$ */
    static double f4_0[2];
    static double f6[2];     /* f_6 */
    static double f7[2];     /* $f_7(step)$ where $step=stepP$ or $stepK$ */
     /* $f_8(\delta,step)$ where $\delta=\delta P_0$ or
      * $\delta k_0$ and $step=stepP$ or $stepK$
      */
    static double f8[2];
    static double fs8[2];    /* $f_8(step,step)$ */
    /* $f_9(\delta)$ where $\delta=\delta p_0$ or $\delta k_0$ */
    static double f9[2];
    static double fs9[2];    /* $f_9(step) $ */
    /* $f1(\delta,m,q)$ where $\delta=\delta k_0$ or $\delta p_0$ */
    static double f1[mqMax2][mqMax4];
    /* $f_2(\delta p,\delta k,q)$ or $f_2(\delta k,\delta p,m$ */
    static double f2[mqMax4];
    /* initial $f_2(\delta p_0,\delta k,q)$ or $f_2(\delta k_0,\delta p,m)$ */
    static double f20[mqMax4];
    static double fs2[mqMax4]; /* $f_2(step,\delta,q)$ */
    static double f3[mqMax4]; /* $f_3(stepK,m)$ or $f_3(stepP,q)$ */
    static double f4[mqMax4]; /* $f_4(\delta k_0,q)$ or $f_4(\delta p_0,m)$ */
/*
 * local variables for the Gaussian
 */
    /*
     * $e^{-\pi\frac{(\delta p+mN)^2}{A}}$ or
     * $e^{-\pi\frac{(\delta k+qN)^2}{B}}$
     */
    static double vg[mqMax2];
    /*
     * $e^{-\pi\frac{(\delta p+mN)^2}{A}}$ or
     * $e^{-\pi\frac{(\delta k+qN)^2}{B}}$
     */
    static double gauss, gx;
/*
 * local variables for parameters
 */
    long N=1<<L;
    long dk, dp;
    int m, q;
    long k2=0, p2=0, k2N=0, p2N=0;
    long k20=0, p20=0, k2N0=0, p2N0=0;
    long pMaxIndx, kMaxIndx;
    long kpIndxReal, kpIndx0Real, kpNIndxReal, kpNIndx0Real;
    long kpIndxImag, kpIndx0Imag, kpNIndxImag, kpNIndx0Imag;
    long kNpIndxReal, kNpIndx0Real, kNpNIndxReal, kNpNIndx0Real;
    long kNpIndxImag, kNpIndx0Imag, kNpNIndxImag, kNpNIndx0Imag;
    long n, minjQ, minjM;
    long mLB, mRB; /* the left bound and right bound for m */
    long qLB, qRB; /* the left bound and the right bound for q */
    double ggR, ggI; /* the real and imaginary parts of the inner product */
    double ggRN, ggIN;
    short dkFirst; /* looping dk first */
    short mFirst;  /* looping m first */
    short dkZero; /* $\delta K_0$ is $0$ or $stepK/2$ */
    short dpZero; /* $\delta P_0$ is $0$ or $stepP/2$ */
    double *vTrans;
    long Ljl2=L-j2+l2, kIndxStep;
/*
 * temporary variables
 */
    double cpi; /* $2\pi 2^{2j_1}/A$ */
    long posIndxR, posIndxI, negIndxR, negIndxI;
    long posIndxR1,posIndxI1,posPrevR,negIndxR1,negIndxI1,negPrevR;
    long mIndx, qIndxR, qIndxI, gIndx, gIndxStep;
    double dtmp, dtmp1, gtmp;
    long ntmp, ntmp1;
/* debug variables */
#ifdef DEBUG0
double dkQN, dpMN;
double gM, gQ;
double dbggR, dbggI;
double dbcR, dbcI;
#endif
/*
 * functions used
 */
    void genericRecursion();
    double gbGaussian();
    long P2(), K2();
/* extracting information from the flag */
    dkFirst = (flag&K_FLAG)==K_FLAG;
    mFirst = (flag&M_FLAG)==M_FLAG;
    dkZero = (flag&A_FLAG)==A_FLAG;
    dpZero = (flag&B_FLAG)==B_FLAG;
/*
 * compute the initial values
 */
    cpi = 2.0*M_PI*(double)(1<<(j1<<1))/(double)A;
    /* compute f4_0 */
    dtmp = cpi*(double)dk0;
    f4_0[0] = cos(dtmp);
    f4_0[1] = -sin(dtmp);
    /* compute f6 */
    dtmp = cpi*(double)N;
    f6[0] = cos(dtmp);
    f6[1] = -sin(dtmp);
    /* compute f7 */
    dtmp = cpi*(double)stepK;
    f7[0] = cos(dtmp);
    f7[1] = -sin(dtmp);
    /* compute f8 */
    dtmp = cpi*(double)dp0*(double)stepK/(double)N;
    f8[0] = cos(dtmp);
    f8[1] = -sin(dtmp);
    dtmp = cpi*(double)stepK*(double)stepP/(double)N;
    fs8[0] = cos(dtmp);
    fs8[1] = -sin(dtmp);
    /* compute f9 */
    dtmp = cpi*dp0;
    f9[0] = cos(dtmp);
    f9[1] = -sin(dtmp);
    dtmp = cpi*(double)stepP;
    fs9[0] = cos(dtmp);
    fs9[1] = -sin(dtmp);
    n = abs(j2-j1);
    pMaxIndx = P2MAX(L,j2,l2);
    kMaxIndx = K2MAX(j2,h2);
    vTrans = trans[j2 - MinOctave + 1]->values;
    kIndxStep = 1<<(Ljl2+1);

/* DEBUG */
#ifdef DEBUG0
fprintf( foutput, "------\nj2 = %d\tl2 = %d\th2 = %d\n",j2,l2,h2);
fprintf( foutput, "dkFirst = %d\tmFirst = %d\n",dkFirst,mFirst);
fprintf( foutput, "L = %i\n",L);
fprintf( foutput, "Ljl2 = %i\n",Ljl2);
if (dkFirst)
{
fprintf( foutput, "stepK = %d\tstepP = %d\n",stepK,stepP);
fprintf( foutput, "dkMaxP = %d\tdpMaxP = %d\n",dkMaxP,dpMaxP);
}
else
{
fprintf( foutput, "stepK = %d\tstepP = %d\n",stepP,stepK);
fprintf( foutput, "dkMaxP = %d\tdpMaxp = %d\n",dpMaxP,dkMaxP);
}
#endif
#ifdef DEBUG0
fprintf( foutput, "f1_0R = %le\tf1_0I = %le\n",f1_0[0],f1_0[1]);
fprintf( foutput, "f2_0R = %le\tf2_0I = %le\n",f2_0[0],f2_0[1]);
fprintf( foutput, "fs2_0R = %le\tfs2_0I = %le\n",fs2_0[0],fs2_0[1]);
fprintf( foutput, "f3_0R = %le\tf3_0I = %le\n",f3_0[0],f3_0[1]);
fprintf( foutput, "f4R = %le\tf4_0I = %le\n",f4[0],f4[1]);
fprintf( foutput, "f6R = %le\tf6_0I = %le\n",f6[0],f6[1]);
fprintf( foutput, "f7R = %le\tf7_0I = %le\n",f7[0],f7[1]);
fprintf( foutput, "f8R = %le\tf8_0I = %le\n",f8[0],f8[1]);
fprintf( foutput, "fs8R = %le\tfs8_0I = %le\n",fs8[0],fs8[1]);
fprintf( foutput, "f9R = %le\tf9_0I = %le\n",f9[0],f9[1]);
fprintf( foutput, "fs9R = %le\tfs9_0I = %le\n",fs9[0],fs9[1]);
#endif
/*
 * Procedure 1: calculation of $f_1(\delta_0,m,q)$
 */
    posIndxR = POSINDREAL; /* q for the real part + */
    negIndxR = NEGINDREAL; /* q for the real part - */
    posIndxI = POSINDIMAG; /* q for the imaginary part + */
    negIndxI = NEGINDIMAG; /* q for the imaginary part - */
    /*
     * initialize the zero position for $f_4$
     */
    f4[mqMax] = f4_0[0];
    f4[mqMax3] = f4_0[1];
    /*
     * initialize the zero position for $f_3$
     */
    f3[mqMax] = f3_0[0];
    f3[mqMax3] = f3_0[1];
    /*
     * initialize the zero position for $f_2$
     */
    f2[mqMax] = f20[mqMax] = f2_0[0];
    f2[mqMax3] = f20[mqMax3] = f2_0[1];
    /*
     * initialize the zero position for $fs_2$
     */
    fs2[mqMax] = fs2_0[0];
    fs2[mqMax3] = fs2_0[1];
    /*
     * initialize the zero position for $f_1$
     */
    /* $f_1($\delta k_0,0,0)$ */
    f1[mqMax][mqMax] = f1_0[0];
    f1[mqMax][mqMax3] = f1_0[1];
    /* $f_1($\delta k_0,m,0)$ */
    posIndxR1 = POSINDREAL; /* m for the real part + */
    negIndxR1 = NEGINDREAL; /* m for the real part - */
    for (m=1;m<mqMax;m++)
	{
	posPrevR = posIndxR1-1;
	negPrevR = negIndxR1+1;
	f1[posIndxR1][mqMax] = CMREAL(f1[posPrevR][mqMax],
		f1[posPrevR][mqMax3],f4[mqMax],f4[mqMax3]);
	f1[posIndxR1][mqMax3] = CMIMAGINARY(f1[posPrevR][mqMax],
		f1[posPrevR][mqMax3],f4[mqMax],f4[mqMax3]);
	f1[negIndxR1][mqMax] = CMREAL(f1[negPrevR][mqMax],
		f1[negPrevR][mqMax3],f4[mqMax],-f4[mqMax3]);
	f1[negIndxR1][mqMax3] = CMIMAGINARY(f1[negPrevR][mqMax],
		f1[negPrevR][mqMax3],f4[mqMax],-f4[mqMax3]);
	posIndxR1++;
	negIndxR1--;
	}
    for (q=1;q<mqMax;q++)
	{
	/*
	 * compute $f_1(\delta k,0,q)
	 */
	f1[mqMax][posIndxR] = f1_0[0];
	f1[mqMax][negIndxR] = f1_0[0];
	f1[mqMax][posIndxI] = f1_0[1];
	f1[mqMax][negIndxI] = f1_0[1];
	/*
	 * calculate $f_4$
	 */
	genericRecursion(f4,posIndxR,posIndxI,negIndxR,negIndxI,f6[0],f6[1]);
	/*
	 * Procedure 2: calculation of $f_3(stepK,m)$ or $f_3(stepP,q)$
	 */
	genericRecursion(f3,posIndxR,posIndxI,negIndxR,negIndxI,f7[0],f7[1]);
	/*
	 * Procedure 4: calculation of $f_2(\delta p_0,\delta k_0,q)$
	 * 		or $f_2(\delta k_0,\delta p_0,m)$
	 */
	genericRecursion(f2,posIndxR,posIndxI,negIndxR,negIndxI,f9[0],f9[1]);
        f20[posIndxR] = f2[posIndxR];
        f20[posIndxI] = f2[posIndxI];
        f20[negIndxR] = f2[negIndxR];
        f20[negIndxI] = f2[negIndxI];
	/*
	 * Procedure 4': calculation of $f_2(stepP,\delta k_0,q)$
	 *		or $f_2(stepK,\delta p_0,m)$
	 */
	genericRecursion(fs2,posIndxR,posIndxI,negIndxR,negIndxI,fs9[0],fs9[1]);
	/*
	 * indeces for looping over m
	 */
	posIndxR1 = POSINDREAL; /* m for the real part + */
	negIndxR1 = NEGINDREAL; /* m for the real part - */

	for (m=1;m<mqMax;m++)
	    {
	    posPrevR = posIndxR1-1;
	    negPrevR = negIndxR1+1;
	    f1[posIndxR1][posIndxR] = CMREAL(f1[posPrevR][posIndxR],
			f1[posPrevR][posIndxI],f4[posIndxR],f4[posIndxI]);
	    f1[posIndxR1][negIndxR] = CMREAL(f1[posPrevR][negIndxR],
			f1[posPrevR][negIndxI],f4[negIndxR],f4[negIndxI]);
	    f1[negIndxR1][negIndxR] = CMREAL(f1[negPrevR][negIndxR],
			f1[negPrevR][negIndxI],f4[negIndxR],-f4[negIndxI]);
	    f1[negIndxR1][posIndxR] = CMREAL(f1[negPrevR][posIndxR],
			f1[negPrevR][posIndxI],f4[posIndxR],-f4[posIndxI]);
	    f1[posIndxR1][posIndxI] = CMIMAGINARY(f1[posPrevR][posIndxR],
			f1[posPrevR][posIndxI],f4[posIndxR],f4[posIndxI]);
	    f1[posIndxR1][negIndxI] = CMIMAGINARY(f1[posPrevR][negIndxR],
			f1[posPrevR][negIndxI],f4[negIndxR],f4[negIndxI]);
	    f1[negIndxR1][negIndxI] = CMIMAGINARY(f1[negPrevR][negIndxR],
			f1[negPrevR][negIndxI],f4[negIndxR],-f4[negIndxI]);
	    f1[negIndxR1][posIndxI] = CMIMAGINARY(f1[negPrevR][posIndxR],
			f1[negPrevR][posIndxI],f4[posIndxR],-f4[posIndxI]);
	    posIndxR1++;
	    negIndxR1--;
	    } /* end of m loop */
	posIndxR++;
	negIndxR--;
	posIndxI++;
	negIndxI--;
	} /* end of q loop */
    /*
     * calculate the initial k2 and p2
     */
    if (dkFirst) {
        minjQ = MIN(L-j1,L-j2);
	//fprintf( foutput, "minjQ = %li\n",minjQ);
        minjM = MIN(j1,j2);
	/* the case of looping $\delta k$ first */
        if (dkZero) {
           ntmp = L-j2-h2;
	   k20 = k2 = K2(dk0,k1,ntmp,N);
           k2N0 = k2N = K2(-dk0,k1,ntmp,N);
        } else {
	    k20 = k2 = K2(dk0,k1,L-j2-h2,N);
        } /* endif */
        if (dpZero) {
           ntmp = j2-l2;
	   p20 = p2 = P2(dp0,p1,ntmp,N);
	   if (p2<0)
		p20 = p2 = pMaxIndx+p2;
           p2N0 = p2N = P2(-dp0,p1,ntmp,N);
	   if (p2N<0)
		p2N0 = p2N = pMaxIndx+p2N;
        } else {
	   p20 = p2 = P2(dp0,p1,j2-l2,N);
	   if (p2<0)
		p20 = p2 = p2+pMaxIndx;
        } /* endif */
    } else {
        minjQ = MIN(j1,j2);
        minjM = MIN(L-j1,L-j2);
	/* the case of looping $\delta p$ first */
        if (dkZero) {
            ntmp = L-j2-h2;
	    k20 = k2 = K2(dp0,k1,ntmp,N);
            k2N0 = k2N = K2(-dp0,k1,ntmp,N);
        } else {
	    k20 = k2 = K2(dp0,k1,L-j2-h2,N);
        } /* endif */
        if (dpZero) {
           ntmp = j2-l2;
	   p20 = p2 = P2(dk0,p1,ntmp,N);
	   if (p2<0)
		p20 = p2 = pMaxIndx+p2;
           p2N0 = p2N = P2(-dk0,p1,ntmp,N);
	   if (p2N<0)
		p2N0 = p2N = pMaxIndx+p2N;
        } else {
	   p20 = p2 = P2(dk0,p1,j2-l2,N);
	   if (p2<0)
		p20 = p2 = pMaxIndx+p2;
        } /* endif */
    } /* end if */
    /* 
     * calculate the indeces for the transformation matrix
     */
    if (k2<0) {
        kpIndxReal = kpIndx0Real = TINDXREAL(-k2,p2,Ljl2);
        kpIndxImag = kpIndx0Imag = TINDXIMAG(-k2,p2,Ljl2);
        if (dpZero) {
            kpNIndxReal = kpNIndx0Real = TINDXREAL(-k2,p2N,Ljl2);
            kpNIndxImag = kpNIndx0Imag = TINDXIMAG(-k2,p2N,Ljl2);
        }
    } else {
        kpIndxReal = kpIndx0Real = TINDXREAL(k2,p2,Ljl2);
        kpIndxImag = kpIndx0Imag = TINDXIMAG(k2,p2,Ljl2);
        if (dpZero) {
            kpNIndxReal = kpNIndx0Real = TINDXREAL(k2,p2N,Ljl2);
            kpNIndxImag = kpNIndx0Imag = TINDXIMAG(k2,p2N,Ljl2);
        }
    } /* endif */
    if (dkZero) {
        if (k2N<0) {
                kNpIndxReal = kNpIndx0Real = TINDXREAL(-k2N,p2,Ljl2);
                kNpIndxImag = kNpIndx0Imag = TINDXIMAG(-k2N,p2,Ljl2);
                if (dpZero) {
                     kNpNIndxReal = kNpNIndx0Real = TINDXREAL(-k2N,p2N,Ljl2);
                     kNpNIndxImag = kNpNIndx0Imag = TINDXIMAG(-k2N,p2N,Ljl2);
                }
        } else {
                kNpIndxReal = kNpIndx0Real = TINDXREAL(k2N,p2,Ljl2);
                kNpIndxImag = kNpIndx0Imag = TINDXIMAG(k2N,p2,Ljl2);
                if (dpZero) {
                    kNpNIndxReal = kNpNIndx0Real = TINDXREAL(k2N,p2N,Ljl2);
                    kNpNIndxImag = kNpNIndx0Imag = TINDXIMAG(k2N,p2N,Ljl2);
                }
        } /* endif */
    }
    /*******************************
     * the generic $\delta k$ loop *
     *******************************/
    for (dk=dk0;dk<=dkMaxP;dk+=stepK)
        {
/* DEBUG */
#ifdef DEBUG0
if (dkFirst)
{
fprintf( foutput, "dk = %d\t",dk);
dtmp = (double)(p1*dk)*2.0*M_PI/(double)N;
dbcR = cos(dtmp);
dbcI = -sin(dtmp);
dtmp = CMREAL(dbcR,dbcI,cR,cI);
dbcI = CMIMAGINARY(dbcR,dbcI,cR,cI);
dbcR = dtmp;
}
else
fprintf( foutput, "dp = %d\t",dk);
#endif
        /*
         * calculate the Gaussian for generic $\delta k$
         */
        ntmp = (-mqMax+1)*N+dk;
        ntmp1 = pnIep[n+1];
        mIndx = 1;
        for (m=1; m<mqMax2; m++) {
           gIndx = GINDX(abs(ntmp),maxlh,minjQ,pnIep[n]);
            if (gIndx>=ntmp1) {
               vg[mIndx++] = 0.0;
            } else {
               vg[mIndx++] = pfG[gIndx];
            } /* endif */
            ntmp += N;
        } /* endfor */
        /* calculate the left and the right bound for q */
        qLB = LB(L,minjQ,pfB[n],dk);
        qRB = RB(L,minjQ,pfB[n],dk);
        /*******************************
         * the generic $\delta p$ loop *
         *******************************/
        for (dp=dp0;dp<=dpMaxP;dp+=stepP) {
            /* calculate the left and the right bound of m */
            mLB = LB(L,minjM,pfB[n],dp);
            mRB = RB(L,minjM,pfB[n],dp);
/* DEBUG */
#ifdef DEBUG0
if (dkFirst)
fprintf( foutput, "dp = %d\n",dp);
else
{
fprintf( foutput, "dk = %d\n",dp);
dtmp = (double)(p1*dp)*2.0*M_PI/(double)N;
dbcR = cos(dtmp);
dbcI = -sin(dtmp);
dtmp = CMREAL(dbcR,dbcI,cR,cI);
dbcI = CMIMAGINARY(dbcR,dbcI,cR,cI);
dbcR = dtmp;
}
fprintf( foutput, "qLB = %d\tqRB = %d\n",qLB,qRB);
fprintf( foutput, "mLB = %d\tmRB = %d\n",mLB,mRB);
#endif
	    /*
             * ---- the double loops m and q ---- 
	     */
            /*
	     * ---- compute the inner product ----
	     */
            if (mFirst)
            {
               /*
		* the case of looping m first
		*/
               ggR = ggI = 0.0;
               /* the case of looping m first */
               mIndx = mqMax+mLB;
               ntmp = dp+mLB*N;
#ifdef DEBUG0
dpMN = (double)ntmp;
dbggR = dbggI = 0.0;
#endif
               ntmp1 = pnIep[n+1];
               for (m=mLB;m<=mRB ; m++) {
                  gIndx = GINDX(abs(ntmp),maxlh,minjM,pnIep[n]);
                  ntmp += N;
                  if (gIndx<ntmp1) {
                     gauss = pfG[gIndx];
                  } else {
#ifdef DEBUG0
gauss = 0.0;
#else
                     mIndx++;
                     continue;
#endif
                  } /* endif */
#ifdef DEBUG0
dkQN = dk+qLB*N;
#endif
                  qIndxR = mqMax+qLB;
                  qIndxI = mqMax3+qLB;
                  for (q=qLB; q<=qRB; q++) {
                     dtmp = gauss*vg[qIndxR];
                     /*
		      * calculate the complex exponential add to
		      * the inner product
		      */
                     ggR += CMREAL(f1[mIndx][qIndxR],f1[mIndx][qIndxI],
				f2[qIndxR],f2[qIndxI])*dtmp;
                     ggI += CMIMAGINARY(f1[mIndx][qIndxR],f1[mIndx][qIndxI],
                                f2[qIndxR],f2[qIndxI])*dtmp;
		     qIndxI++;
		     qIndxR++;
#ifdef DEBUG0
dtmp = cpi*dpMN*dkQN/(double)N;
if (dkFirst)
dtmp1 = exp(-M_PI*((double)(dkQN*dkQN)/(double)B+
	(double)(dpMN*dpMN)/(double)A));
else
dtmp1 = exp(-M_PI*((double)(dkQN*dkQN)/(double)A+
	(double)(dpMN*dpMN)/(double)B));
dbggR += cos(dtmp)*dtmp1;
dbggI -= sin(dtmp)*dtmp1;
dkQN += (double)N;
#endif
                  } /* endfor q */
		  mIndx++;
#ifdef DEBUG0
dpMN += (double)N;
#endif
               } /* endfor m */
            } else {
	       /*
		* the case of looping q first
		*/
               ggR = ggI = 0.0;
               /* the case of looping m first */
               qIndxR = mqMax+qLB;
               qIndxI = mqMax3+qLB;
#ifdef DEBUG0
dkQN = (double)(dk+qLB*N);
dbggR = dbggI = 0.0;
#endif
               for (q=qLB;q<=qRB ; q++) {
                  ntmp = dp+mLB*N;
                  ntmp1 = pnIep[n+1];
                  gtmp = vg[qIndxR];
                  mIndx = mqMax+mLB;
#ifdef DEBUG0
dpMN = (double)ntmp;
#endif
                  for (m=mLB; m<=mRB; m++) {
                    gIndx = GINDX(abs(ntmp),maxlh,minjM,pnIep[n]);
                     /* calculate the gassian */
                     ntmp+=N;
                     if (gIndx<ntmp1) {
                        gauss = pfG[gIndx]*gtmp;
                     } else {
#ifdef DEBUG0
gauss = 0.0;
#else
                        mIndx++;
                        continue;
#endif
                     } /* endif */
                     /*
		      * calculate the complex exponential add
		      * to the inner product
		      */
                     ggR += CMREAL(f1[mIndx][qIndxR],f1[mIndx][qIndxI],
				f2[qIndxR],f2[qIndxI])*gauss;
                     ggI += CMIMAGINARY(f1[mIndx][qIndxR],f1[mIndx][qIndxI],
                                f2[qIndxR],f2[qIndxI])*gauss;
		     mIndx++;
#ifdef DEBUG0
dtmp = cpi*dpMN*dkQN/(double)N;
if (dkFirst)
dtmp1 = exp(-M_PI*((double)(dkQN*dkQN)/(double)B+
                (double)(dpMN*dpMN)/(double)A));
else
dtmp1 = exp(-M_PI*((double)(dkQN*dkQN)/(double)A+
                (double)(dpMN*dpMN)/(double)B));
dbggR += cos(dtmp)*dtmp1;
dbggI -= sin(dtmp)*dtmp1;
dpMN += (double)N;
#endif
                  } /* endfor m */
		  qIndxR++;
		  qIndxI++;
#ifdef DEBUG0
dkQN += (double)N;
#endif
               } /* endfor q */
            } /* endif */
	/************************************
	 * update the transfromation matrix *
	 ************************************/
	if (k2<0)
	    {
	    vTrans[kpIndxReal] -= ggR;
	    vTrans[kpIndxImag] += ggI;
	    }
        else if (k2==0 || k2==kMaxIndx)
            {
             vTrans[kpIndxReal] -= ggR+ggR;
             }
	else
	    {
	    vTrans[kpIndxReal] -= ggR;
	    vTrans[kpIndxImag] -= ggI;
	    }
	if (dkZero||dpZero)
	    {
	    ggRN = CMREAL(ggR,-ggI,cos_2phi,sin_2phi);
	    ggIN = CMIMAGINARY(ggR,-ggI,cos_2phi,sin_2phi);
	    }
	if (dkZero&&k2N!=k2&&k2N!=-kMaxIndx)
	    {
	    if (k2N<0)
		{
	    	vTrans[kNpIndxReal] -= ggRN;
	    	vTrans[kNpIndxImag] += ggIN;
		}
             else if (k2N==0 || k2N==kMaxIndx)
                {
                 vTrans[kNpIndxReal] -= ggRN+ggRN;
                }
	    else
	    	{
	    	vTrans[kNpIndxReal] -= ggRN;
	    	vTrans[kNpIndxImag] -= ggIN;
	    	}
	    }
	if (dpZero&&p2N!=p2)
	    {
	    if (dkFirst)
		ntmp = dk*p1*2;
	    else
		ntmp = 2*dp*p1;
	    mIndx = abs(ntmp)%N;
	    dtmp = pfCE1[mIndx];
	    if (ntmp>0)
		dtmp1 = -pfCE1[N+mIndx];
	    else
		dtmp1 = pfCE1[N+mIndx];
	    if (k2<0)
		{
		vTrans[kpNIndxReal] -= CMREAL(dtmp,dtmp1,ggRN,ggIN);
		vTrans[kpNIndxImag] += CMIMAGINARY(dtmp,dtmp1,ggRN,ggIN);
		}
	    else if (k2 == 0 || k2 == kMaxIndx)
                {
                vTrans[kpNIndxReal] -= 2.0*CMREAL(dtmp,dtmp1,ggRN,ggIN);
                }
            else if (k2!=0)
		{
		vTrans[kpNIndxReal] -= CMREAL(dtmp,dtmp1,ggRN,ggIN);
		vTrans[kpNIndxImag] -= CMIMAGINARY(dtmp,dtmp1,ggRN,ggIN);
		}
	    if (dkZero&&k2N!=k2&&k2N!=-kMaxIndx)
		{
	        if (k2N<0)
		    {
		    vTrans[kNpNIndxReal] -= CMREAL(dtmp,-dtmp1,ggR,ggI);
		    vTrans[kNpNIndxImag] += CMIMAGINARY(dtmp,-dtmp1,ggR,ggI);
		    }
                else if (k2N==0)
                    {
                    vTrans[kNpNIndxReal] -= 2.0*CMREAL(dtmp,-dtmp1,ggR,ggI);
                    }
	    	else
		    {
		    vTrans[kNpNIndxReal] -= CMREAL(dtmp,-dtmp1,ggR,ggI);
		    vTrans[kNpNIndxImag] -= CMIMAGINARY(dtmp,-dtmp1,ggR,ggI);
		    }
		}
	    }
/* DEBUG */
#ifdef DEBUG0
fprintf( foutput, "k2 = %d\tp2 = %d\n",k2,p2);
if (dpZero&&p2N!=p2)
      fprintf( foutput, "k2 = %d\tp2N = %d\n",k2,p2N);
if (dkZero&&k2N!=k2) {
      fprintf( foutput, "k2N = %d\tp2 = %d\n",k2N,p2);
   if (dpZero&&p2N!=p2)
      fprintf( foutput, "k2N = %d\tp2N = %d\n",k2N,p2N);
}
fprintf( foutput, "kpIndxReal = %d\tkpIndxImag = %d\n",kpIndxReal,kpIndxImag);
if (dpZero&&p2N!=p2)
      fprintf( foutput, "kpNIndxReal = %d\tkpNIndxImag = %d\n",kpNIndxReal,kpNIndxImag);
if (dkZero&&k2N!=k2) {
      fprintf( foutput, "kNpIndxReal = %d\tkNpIndxImag = %d\n",kNpIndxReal,kNpIndxImag);
   if (dpZero&&p2N!=p2)
      fprintf( foutput, "kNpNIndxReal = %d\tkNpNIndxImag = %d\n",kNpNIndxReal,kNpNIndxImag);
}
dtmp = CMREAL(dbggR,dbggI,dbcR,dbcI);
dbggI = CMIMAGINARY(dbggR,dbggI,dbcR,dbcI);
dbggR = dtmp;
fprintf( foutput, "ggR = %lf\tggI = %lf\n",ggR,ggI);
fprintf( foutput, "dbggR = %lf\tdbggI = %lf\n",dbggR,dbggI);
#endif
	if (dp > dpMaxP-stepP)
	    break; /* this is the last step */
	/*
	 * update the necessary variables for the next step
	 */
        /*
         * calculate $k_2$ or $p_2$
         */
         if (!dkFirst) {
	    /* the case of looping $\delta p$ first */
            if (k2==kMaxIndx) {
               k2 = -kMaxIndx+1;
               kpIndxReal -= kIndxStep;
               kpIndxImag -= kIndxStep;
               if (dpZero) {
                   kpNIndxReal -= kIndxStep;
                   kpNIndxImag -= kIndxStep;
               } /* dpZero */
            } else if (k2<0) {
               kpIndxReal -= kIndxStep;
               kpIndxImag -= kIndxStep;
               if (dpZero) {
                      kpNIndxReal -= kIndxStep;
                      kpNIndxImag -= kIndxStep;
               } /* dpZero */
               k2++;
            } else {
               kpIndxReal += kIndxStep;
               kpIndxImag += kIndxStep;
               if (dpZero) {
                      kpNIndxReal += kIndxStep;
                      kpNIndxImag += kIndxStep;
               } /* dpZero */
               k2++;
            } /* k2==kMaxIndx */
            if (dkZero)
                {
		k2N--;
                if (k2N<0) {
                   kNpIndxReal += kIndxStep;
                   kNpIndxImag += kIndxStep;
                   if (dpZero) {
                      kNpNIndxReal += kIndxStep;
                      kNpNIndxImag += kIndxStep;
                   } /* dpZero */
                } else {
                   kNpIndxReal -= kIndxStep;
                   kNpIndxImag -= kIndxStep;
                   if (dpZero) {
                      kNpNIndxReal -= kIndxStep;
                      kNpNIndxImag -= kIndxStep;
                   } /* dpZero */
                } /* k2N<0 */
                } /* dkZero */
         } else {
	    /* the case of looping $\delta k$ first */
            if (p2==pMaxIndx-1) {
               p2 = 0;
               kpIndxReal -= pMaxIndx-1;
               kpIndxImag -= pMaxIndx-1;
               if (dkZero) {
                  kNpIndxReal -= pMaxIndx-1;
                  kNpIndxImag -= pMaxIndx-1;
               } /* dkZero */
            } else {
               p2++;
               kpIndxReal++;
               kpIndxImag++;
               if (dkZero) {
                  kNpIndxReal++;
                  kNpIndxImag++;
               } /* dkZero */
            } /* p2==pMaxIndx-1 */
            if (dpZero) {
                  if (p2N==0) {
                        p2N = pMaxIndx-1;
                        kpNIndxReal += p2N;
                        kpNIndxImag += p2N;
                        if (dkZero) {
                            kNpNIndxReal += p2N;
                            kNpNIndxImag += p2N;
                        } /* dkZero */
              } else {
                        p2N--;
                        kpNIndxReal--;
                        kpNIndxImag--;
                        if (dkZero) {
                            kNpNIndxReal--;
                            kNpNIndxImag--;
                        } /* dkZero */
                   } /* p2N==0 */
              } /* dpZero */
         } /* !dkFirst */
	 /*
	  * Procedure 6: calculation of $f_2(\delta p,\delta k, q)$
	  */
	 posIndxR = POSINDREAL;
	 posIndxI = POSINDIMAG;
	 negIndxR = NEGINDREAL;
	 negIndxI = NEGINDIMAG;
         dtmp = CMREAL(f2[POSZEROREAL],f2[POSZEROIMAG],
                       fs2[POSZEROREAL],fs2[POSZEROIMAG]);
         f2[POSZEROIMAG] = CMIMAGINARY(f2[POSZEROREAL],f2[POSZEROIMAG],
                       fs2[POSZEROREAL],fs2[POSZEROIMAG]);
         f2[POSZEROREAL] = dtmp;
	 for (q=1;q<mqMax;q++)
	     gnrRcsArray(f2,
		   posIndxR++,posIndxI++,negIndxR--,negIndxI--,
		   fs2);
	} /* end of dp loop */
	if (dk > dkMaxP-stepK)
	    break; /* this is the last step */
	/*
	 * update the variables for the next step
	 */
        /*
         * calculate $k_2$ or $p_2$
         */
        if (dkFirst) {
	   /* the case of looping $\delta k$ first */
           p2 = p20;
           p2N = p2N0;
            if (k2==kMaxIndx) {
               k2 = -kMaxIndx+1;
               kpIndxReal = kpIndx0Real = kpIndx0Real-kIndxStep;
               kpIndxImag = kpIndx0Imag = kpIndx0Imag-kIndxStep;
               if (dpZero) {
                   kpNIndxReal = kpNIndx0Real = kpNIndx0Real-kIndxStep;
                   kpNIndxImag = kpNIndx0Imag = kpNIndx0Imag-kIndxStep;
               } /* dpZero */
            } else if (k2<0) {
               kpIndxReal = kpIndx0Real = kpIndx0Real-kIndxStep;
               kpIndxImag = kpIndx0Imag = kpIndx0Imag-kIndxStep;
               if (dpZero) {
                      kpNIndxReal = kpNIndx0Real = kpNIndx0Real-kIndxStep;
                      kpNIndxImag = kpNIndx0Imag = kpNIndx0Imag-kIndxStep;
               } /* dpZero */
               k2++;
            } else {
               kpIndxReal = kpIndx0Real = kpIndx0Real+kIndxStep;
               kpIndxImag = kpIndx0Imag = kpIndx0Imag+kIndxStep;
               if (dpZero) {
                      kpNIndxReal = kpNIndx0Real = kpNIndx0Real+kIndxStep;
                      kpNIndxImag = kpNIndx0Imag = kpNIndx0Imag+kIndxStep;
               } /* dpZero */
               k2++;
            } /* k2==kMaxIndx */
            if (dkZero)
                {
		k2N--;
                if (k2N<0) {
                   kNpIndxReal = kNpIndx0Real = kNpIndx0Real+kIndxStep;
                   kNpIndxImag = kNpIndx0Imag = kNpIndx0Imag+kIndxStep;
                   if (dpZero) {
                      kNpNIndxReal = kNpNIndx0Real = kNpNIndx0Real+kIndxStep;
                      kNpNIndxImag = kNpNIndx0Imag = kNpNIndx0Imag+kIndxStep;
                   } /* dpZero */
                } else {
                   kNpIndxReal = kNpIndx0Real = kNpIndx0Real-kIndxStep;
                   kNpIndxImag = kNpIndx0Imag = kNpIndx0Imag-kIndxStep;
                   if (dpZero) {
                      kNpNIndxReal = kNpNIndx0Real = kNpNIndx0Real-kIndxStep;
                      kNpNIndxImag = kNpNIndx0Imag = kNpNIndx0Imag-kIndxStep;
                   } /* dpZero */
                } /* k2N<0 */
                } /* dkZero */
        } else { 
	   /* the case of looping $\delta p$ first */
           k2 = k20;
           k2N = k2N0;
            if (p2==pMaxIndx-1) {
               p2 = 0;
               kpIndxReal = kpIndx0Real = kpIndx0Real-pMaxIndx+1;
               kpIndxImag = kpIndx0Imag = kpIndx0Imag-pMaxIndx+1;
               if (dkZero) {
                  kNpIndxReal = kNpIndx0Real = kNpIndx0Real-pMaxIndx+1;
                  kNpIndxImag = kNpIndx0Imag = kNpIndx0Imag-pMaxIndx+1;
               } /* dkZero */
            } else {
               p2++;
               kpIndxReal = ++kpIndx0Real;
               kpIndxImag = ++kpIndx0Imag;
               if (dkZero) {
                  kNpIndxReal = ++kNpIndx0Real;
                  kNpIndxImag = ++kNpIndx0Imag;
               } /* dkZero */
            } /* p2==pMaxIndx-1 */
            if (dpZero) {
                  if (p2N==0) {
                        p2N = pMaxIndx-1;
                        kpNIndxReal = kpNIndx0Real = kpNIndx0Real+p2N;
                        kpNIndxImag = kpNIndx0Imag = kpNIndx0Imag+p2N;
                        if (dkZero) {
                            kNpNIndxReal = kNpNIndx0Real = kNpNIndx0Real+p2N;
                            kNpNIndxImag = kNpNIndx0Imag = kNpNIndx0Imag+p2N;
                        } /* dkZero */
                  } else {
                        p2N--;
                        kpNIndxReal = --kpNIndx0Real;
                        kpNIndxImag = --kpNIndx0Imag;
                        if (dkZero) {
                            kNpNIndxReal = --kNpNIndx0Real;
                            kNpNIndxImag = --kNpNIndx0Imag;
                        } /* dkZero */
                   } /* p2N==0 */
              } /* dpZero */
        } /* endif */
	/*
	 * Procedure 3: calculation of $f_1(\delta k,m,q)$
	 */
        /* update $f_1(\delta k,0,0) */
        dtmp = CMREAL(f1[POSZEROREAL][POSZEROREAL],f1[POSZEROREAL][POSZEROIMAG],
                                        f3[POSZEROREAL],f3[POSZEROIMAG]);
        f1[POSZEROREAL][POSZEROIMAG] = CMIMAGINARY(
                      f1[POSZEROREAL][POSZEROREAL],f1[POSZEROREAL][POSZEROIMAG],
                      f3[POSZEROREAL],f3[POSZEROIMAG]);
        f1[POSZEROREAL][POSZEROREAL] = dtmp;
        /* update $f_1(\delta k,0,q) */
	posIndxR = POSINDREAL; /* index for positive q of the real part */
	posIndxI = POSINDIMAG;    /* index for positive q of the imaginary part */
	negIndxR = NEGINDREAL; /* index for negative q of the real part */
	negIndxI = NEGINDIMAG;   /* index for negative q of the imaginary part */
        for (q=1; q<mqMax; q++) {
           gnrRcs(f1[POSZEROREAL],
                  posIndxR++,posIndxI++,negIndxR--,negIndxI--,f3[POSZEROREAL],
                  f3[POSZEROIMAG]);
        } /* endfor */
        /* update $f_1(\delta k,m,q)$ for $m>1$ */
        posIndxR1 = POSINDREAL;  /* index for positive m of the real part */
        posIndxI1 = POSINDIMAG;    /* index for positive m of the imaginary part */
        negIndxR1 = NEGINDREAL; /* index for negative m of the real part */
        negIndxI1 = NEGINDIMAG;    /* index for negative m of the imaginary part */
	for (m=1;m<mqMax;m++)
            {
	    posIndxR = POSINDREAL; /* index for positive q of the real part */
	    posIndxI = POSINDIMAG;    /* index for positive q of the imaginary part */
	    negIndxR = NEGINDREAL; /* index for negative q of the real part */
	    negIndxI = NEGINDIMAG;   /* index for negative q of the imaginary part */
            /* update $f_1(\delta k,m,0) */
            dtmp = CMREAL(f1[posIndxR1][POSZEROREAL],f1[posIndxR1][POSZEROIMAG],
                          f3[posIndxR1],f3[posIndxI1]);
            f1[posIndxR1][POSZEROIMAG] = CMIMAGINARY(
                          f1[posIndxR1][POSZEROREAL],f1[posIndxR1][POSZEROIMAG],
                          f3[posIndxR1],f3[posIndxI1]);
            f1[posIndxR1][POSZEROREAL] = dtmp;
            dtmp = CMREAL(f1[negIndxR1][POSZEROREAL],f1[negIndxR1][POSZEROIMAG],
                          f3[negIndxR1],f3[negIndxI1]);
            f1[negIndxR1][POSZEROIMAG] = CMIMAGINARY(
                          f1[negIndxR1][POSZEROREAL],f1[negIndxR1][POSZEROIMAG],
                          f3[negIndxR1],f3[negIndxI1]);
            f1[negIndxR1][POSZEROREAL] = dtmp;
	    for (q=1;q<mqMax;q++)
		{
		gnrRcs(f1[posIndxR1],
			posIndxR,posIndxI,negIndxR,negIndxI,
			f3[posIndxR1],f3[posIndxI1]);
		gnrRcs(f1[negIndxR1],
			posIndxR,posIndxI,negIndxR,negIndxI,
			f3[negIndxR1],f3[negIndxI1]);
		posIndxR++;
		posIndxI++;
		negIndxR--;
		negIndxI--;
		} /* end of q loop */
            posIndxR1++;
            posIndxI1++;
            negIndxR1--;
            negIndxI1--;
            } /* end of m loop */
	/*
	 * Procedure 5 & 5': 	$f_2(\delta p_0,\delta k, q)$ and
	 *			$f_2(stepP,\delta k, q$
	 */
	posIndxR = POSINDREAL;
	posIndxI = POSINDIMAG;
	negIndxR = NEGINDREAL;
	negIndxI = NEGINDIMAG;
        f2[POSZEROREAL] = CMREAL(f20[POSZEROREAL],f20[POSZEROIMAG],
                                        f8[0],f8[1]);
        f2[POSZEROIMAG] = CMIMAGINARY(f20[POSZEROREAL],f20[POSZEROIMAG],
                                        f8[0],f8[1]);
        f20[POSZEROREAL] = f2[POSZEROREAL];
        f20[POSZEROIMAG] = f2[POSZEROIMAG];
        dtmp = CMREAL(fs2[POSZEROREAL],fs2[POSZEROIMAG],
                      fs8[0],fs8[1]);
        fs2[POSZEROIMAG] = CMIMAGINARY(fs2[POSZEROREAL],fs2[POSZEROIMAG],
                      fs8[0],fs8[1]);
        fs2[POSZEROREAL] = dtmp;
	for (q=1;q<mqMax;q++)
	    {
	    gnrRcs(f20,
		posIndxR,posIndxI,negIndxR,negIndxI,
		f8[0],f8[1]);
             f2[posIndxR] = f20[posIndxR];
             f2[posIndxI] = f20[posIndxI];
             f2[negIndxR] = f20[negIndxR];
             f2[negIndxI] = f20[negIndxI];
	    gnrRcs(fs2,
		posIndxR,posIndxI,negIndxR,negIndxI,
		fs8[0],fs8[1]);
	    posIndxR++;
	    posIndxI++;
	    negIndxR--;
	    negIndxI--;
	    } /* end of q loop */
	} /* end of dk loop */
}
/*
 * compute $k_2$ or $p_2$ from delta
 *
 * inputs:
 *	delta:	$\delta k$ or $\delta p$
 *	p1:	$2^{j_1-l_1}p_1$ (long)
 *	j2:		$j_2-l_2$ (long)
 *	N:		signal size (long)
 *
 * output:
 *      p2:          $p_2$ (dobule)
 *
 *      Bugs:
 *              Should be written as macro
*/
long P2(delta,p1,j2,N)
long delta;
long p1;
long j2;
long N;
{
        long tmp;

        tmp = delta+p1;
        if (SIGN(tmp)>0) {
              return((tmp%N)>>j2);
        } else {
              return(-((abs(tmp)%N)>>j2));
        } /* endif */
}
/*
 * compute $k_2$ or $p_2$ from delta
 *
 * inputs:
 *	delta:	$\delta k$ or $\delta p$
 *	k1:	$2^{L-j_1-h_1}k_1$ (long)
 *	Lj2:	$L_2-j_2-h_2$ (long)
 *	N:		signal size (long)
 *
 * output:
 *      k2:          $k_2$ (double)
 *
 *      Bugs:
 *              Should be written as macro
*/
long K2(delta,k1,Lj2,N)
long delta;
long k1;
long Lj2;
long N;
{
        long tmp;

        tmp = delta+k1;
        while (tmp<0) {
           tmp += N;
        } /* endwhile */
        tmp %= N;
        if (tmp>(N>>1)) {
             tmp = N-tmp;
             return(-(tmp>>Lj2));
         } else {
             return(tmp>>Lj2);
         } /* endif */
}
/*
 * compute the initial values of $f_1$, $f_2$, $f_3$
 */
void initValues(f1_0,f2_0,fs2_0,f3_0,dk0,dp0,p1,stepK,stepP,N,
		cpi,pfCE1,cos_phi,sin_phi,flag)
double f1_0[2];
double f2_0[2];
double fs2_0[2];
double f3_0[2];
long dk0, dp0;
long p1;
long stepK, stepP;
long N;
double cpi;
double *pfCE1;
double cos_phi, sin_phi;
unsigned long flag;
{
    long ntmp;
    double dtmp,dtmp1;

    if ((flag&K_FLAG)==K_FLAG)
	   {
	   ntmp = abs(p1*dk0)%N;
           dtmp = pfCE1[ntmp];
           dtmp1 = pfCE1[N+ntmp];
           if (dk0<0) {
              f1_0[0] = CMREAL(dtmp,dtmp1,cos_phi,sin_phi);
              f1_0[1] = CMIMAGINARY(dtmp,dtmp1,cos_phi,sin_phi);
           } else {
              f1_0[0] = CMREAL(dtmp,-dtmp1,cos_phi,sin_phi);
              f1_0[1] = CMIMAGINARY(dtmp,-dtmp1,cos_phi,sin_phi);
           } /* endif */
	   /*
	    * initial values for Procedure 2
	    */
	   /* compute f3_0 */
	   ntmp = (p1*stepK)%N;
	   f3_0[0] = pfCE1[ntmp];
	   f3_0[1] = -pfCE1[ntmp+N];
	   /*
	    * initial values for Procedure 4
	    */
	   /* compute f2_0 */
	   dtmp = cpi*(double)dp0*(double)dk0/(double)N;
	   f2_0[0] = cos(dtmp);
	   f2_0[1] = -sin(dtmp);
	   /*
	    * initial values for Procedure 4'
	    */
	   /* compute f2_0 */
	   dtmp = cpi*(double)dk0*(double)stepP/(double)N;
	   fs2_0[0] = cos(dtmp);
	   fs2_0[1] = -sin(dtmp);
	   }
    else
	   {
	   f1_0[0] = cos_phi;
	   f1_0[1] = sin_phi;
	   /*
	    * initial values for Procedure 2
	    */
	   /* compute f3_0 */
	   f3_0[0] = 1.0;
	   f3_0[1] = 0.0;
	   /*
	    * initial values for Procedure 4
	    */
	   /* compute f2_0 */
	   dtmp = (cpi*(double)dp0+
		2.0*M_PI*(double)p1)*(double)dk0/(double)N;
	   f2_0[0] = cos(dtmp);
	   f2_0[1] = -sin(dtmp);
	   /*
	    * initial values for Procedure 4'
	    */
	   /* compute f2_0 */
	   dtmp = (cpi*(double)dp0+
		2.0*M_PI*(double)p1)*(double)stepK/(double)N;
	   fs2_0[0] = cos(dtmp);
	   fs2_0[1] = -sin(dtmp);
	   }
}
/* the left bound of m or q */
/*
int LB(L,dMax,d)
int L;
long dMax;
long d;
{
   long tmp=dMax+d, tmp1;

   tmp1 = tmp>>(L-1);
   if (tmp1&1) {
      return(-(tmp>>L));
   } else {
      return(-(tmp>>L)-1);
   }
}
*/
/* the right bound of m or q */
/*
int RB(L,dMax,d)
int L;
long dMax;
long d;
{
   long tmp=dMax-d, tmp1;

   tmp1 = tmp >> (L-1);
   if (tmp1&1) {
      return(tmp>>L);
   } else {
      return((tmp>>L)+1);
   }
}
*/
/*
 * end of update.c
 */
