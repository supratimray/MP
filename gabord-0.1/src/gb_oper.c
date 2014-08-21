#include "mpp.h"
/* real part of complex multiplication */
#define CMREAL(a1,b1,a2,b2)     ((a1)*(a2)-(b1)*(b2))
/* imaginary part of complex multiplication */
#define CMIMAGINARY(a1,b1,a2,b2)        ((a1)*(b2)+(a2)*(b1))
/* calculate the index for the transform given by $k_2$ and $p_2$ */
#define TINDXREAL(k2,p2,Ljl2)   (((k2)<<((Ljl2)+1))+(p2))
#define TINDXIMAG(k2,p2,Ljl2)	(((k2)<<((Ljl2)+1))+(1<<(Ljl2))+(p2))
#define SSFACTOR(j,l)   (((j)>(l))?(l):(j)) /* subsample factor */
#define STEP(j,l)       (1<<((j)-(l)))
/* global variables */
extern double *cur_norm;
extern double *pfG, *pfCE1;
extern double *pfCE2, *pfCE3, *pfC;
extern double *pfB;
extern int *pnAep, *pnIep;
extern int cur_l, cur_h; /* defined in ng_cmd.c */
extern double epslonG; /* defined in ng_cmd.c */
extern int nBoundG; /* defined in ng_cmd.c */
INDEX indexG1=(INDEX)NULL, indexG2=(INDEX)NULL;
/*
 * get the maximum modula from a trans
 *
 * Inputs:
 *	trans			transformation from the gabor library
 *
 *	num_octave		number of octaves in the transform
 *
 *	SubsampleOctaveTime	the octave number that begin to subsample 
 *				in translation
 *
 *	SubsampleOctaveFreq	the octave number that begin to subsample
 *				in frequency
 * Bugs:
 *	1. requires the signal size is a power of 2
 *	2. it does not check the validity of 'MacOctave', 'ShiftOctave',
 *	   'SubsampleOctaveTime', 'SubsampleOctaveFreq'.
 */
WRD GaborGetMaxFrmTrans(trans,filter, MinOctave, MaxOctave, num_octave,SubsampleOctaveTime,
		SubsampleOctaveFreq)
SIGNAL *trans;
FILTER *filter;
int MinOctave, MaxOctave; /* Min and Max octave for the decomposition on Gabor functions */
int num_octave;
int SubsampleOctaveTime;
int SubsampleOctaveFreq;
{
    WRD wrd=(WRD)NULL;
    int k, j, SampleRate_f, SampleRate_n;
    long SigSize, HalfSigSize, DoubleSigSize;
    int freq, bndBadK;
    int n_length, n_2_length, kend, index;
    double modula, value_r, value_i; //, nrm;
    double MaxModula = 0.0, MaxOct, MaxFreq, MaxTranslation;
    //double MaxModulaG, MaxOctG, MaxFreqG, MaxTranslationG;
    double MaxReal, MaxImg, MaxNorm;
    //double MaxRealG, MaxImgG, MaxNormG;
    double *pointer;
    int LnSigSize;
    double GaborGetGaborNorm();
    WRD AllocWord();
    void SigAbsMx(), GaborArrayMax();
    extern void complex_array_max();

    if (trans == (SIGNAL *)NULL || filter == (FILTER *)NULL)
	perror("GaborGetMaxFrmTrans(): null argument!");

    if (trans[0] == (SIGNAL)NULL)
	perror("Internal error!");
    DoubleSigSize = trans[0]->size;
    SigSize = DoubleSigSize>>1;
    HalfSigSize = SigSize>>1;
    LnSigSize = (int)log2((double)SigSize);
/*
 * searching maximum modula
 */

/* the Dirac basis: j = 0; index 0 in the array trans */
if ( (MinOctave == 1) && (MaxOctave == LnSigSize - 1) )
	{
	SigAbsMx(trans[0],&MaxReal,&MaxTranslation);
	MaxOct = 0.0;
	MaxFreq = 0.0;
	MaxModula = fabs(MaxReal*MaxReal)/2.0;
	MaxImg = 0.0;
	MaxNorm = 1.0;
	}

/* The Gabor Basis : j = MinOctave, ..., MaxOctave; index in the array trans: j - MinOctave + 1 */
    for (j=MinOctave;j<=MaxOctave;j++)
	{
	if (trans[j - MinOctave + 1] == (SIGNAL)NULL)
	    perror("Internal error!");
	if (j>LnSigSize-SubsampleOctaveFreq)
	    {
	    kend = 1<<(LnSigSize-1);
	    bndBadK = 1<<(LnSigSize-j);
	    SampleRate_f = 1;
	    }
	else
	    {
	    kend = 1<<(SubsampleOctaveFreq+j-2);
	    bndBadK = 1<<(SubsampleOctaveFreq-1);
	    SampleRate_f = 1<<(LnSigSize-SubsampleOctaveFreq-j+1);
	    }
	if (j>=SubsampleOctaveTime)
	    {
	    n_length = SigSize>>(j-SubsampleOctaveTime+1);
	    n_2_length = n_length<<1;
	    SampleRate_n = 1<<(j-SubsampleOctaveTime+1);
	    }
	else
	    {
	    n_length = SigSize;
	    n_2_length = DoubleSigSize;
	    SampleRate_n = 1;
	    }
	pointer = trans[j - MinOctave + 1]->values;
	for (k=0;k<=kend;k++)
	    {
	    if ((k<bndBadK && k>0)||(k>(kend-bndBadK) && k<kend))
		{
		pointer += n_2_length;
		continue;
		}
	    freq = k*SampleRate_f;
	    GaborArrayMax(filter,pointer,n_length,SigSize, LnSigSize,
			j,freq,
			SampleRate_n,
			&modula,&value_r,&value_i,&index);
	    if (freq==0 || freq==HalfSigSize)
		modula /= 2.0;
	    if (modula>MaxModula)
		{
		MaxOct = (double)j;
		MaxFreq = (double)freq;
		MaxTranslation = (double)(SampleRate_n*index);
		MaxModula = modula;
		MaxReal = value_r;
		MaxImg = value_i;
		}
	    pointer += n_2_length;
	    } /* end of loop k */
	} /* end of loop j */


    /* The Fourier Basis : j = LnSigSize; index MaxOctave-MinOctave+2 in the array trans*/
if ( (MinOctave == 1) && (MaxOctave == LnSigSize - 1) )
	{
    	if (trans[MaxOctave-MinOctave+2] == (SIGNAL)NULL)
		perror("Internal error!");
    	pointer = trans[MaxOctave-MinOctave+2]->values;
    	complex_array_max(pointer,SigSize,(SigSize>>1)+1,
		&modula,&value_r,&value_i,&index);
    	if (index==0 || index==HalfSigSize)
		modula /= 2.0;
    	if (modula>MaxModula)
		{
		MaxModula = modula;
		MaxOct = (double)LnSigSize;
		MaxFreq = (double)index;
		MaxTranslation = 0.0;
		MaxReal = value_r;
		MaxImg = value_i;
		}
	}
/*
 * install the maximum information into wrd
 */
    wrd = AllocWord();
    if (MaxOct == (double)LnSigSize) /* Fourier basis */
	{
	if (MaxFreq==0.0 || MaxFreq==(double)HalfSigSize)
	    { 
	    MaxNorm = 1.0/sqrt((double)SigSize);
	    wrd->value = MaxNorm;
	    }
	else
	    {
	    MaxNorm = M_SQRT2/sqrt((double)SigSize);
	    wrd->value = MaxNorm;
	    } 
	}
    else if (MaxOct==0) /* Dirac basis */
        wrd->value = MaxNorm;
    else /* Gabor functions */
	{
	MaxNorm = GaborGetGaborNorm(
			MaxModula,MaxReal,MaxImg,SigSize,
			LnSigSize,(int)MaxOct,(long)MaxFreq,
			(long)MaxTranslation);
	wrd->value = MaxNorm;
	}


    if (MaxOct==0.0 || MaxFreq==0.0 || MaxFreq==(double)HalfSigSize)
	MaxModula = sqrt(2.0*MaxModula);
    else
	{
	if (MaxOct != (double)LnSigSize)
	    {
	    MaxModula = sqrt(MaxModula)*MaxNorm;
	    MaxReal *= MaxNorm;
	    }
	else
	    {
	    MaxModula = M_SQRT2*sqrt(MaxModula);
	    MaxReal *= M_SQRT2;
	    }
	}
    wrd->coeff = MaxModula;
    if (MaxImg>0)
	wrd->index->phase = acos(MaxReal/MaxModula);
    else
	wrd->index->phase = -acos(MaxReal/MaxModula);
    wrd->index->id = MaxFreq;
    wrd->index->octave = MaxOct;
    wrd->index->position = MaxTranslation;

    return(wrd);
}
/*
 *
 * Get residue for trans from a wrd
 *
 * Inputs:
 *	trans		transformation for the gabor library (SIGNAL *)
 *	filter		the basic gabor functions	(FILTER *)
 *	wrd		wrd selected from trans (WRD)
 *
 * Output:
 *	trans[0]	residue of the project persuit (SIGNAL)
 *
 */
void GaborGetResidue(trans,filter,wrd,num_octave)
SIGNAL *trans;
FILTER *filter;
WRD wrd;
int num_octave;
{
    int i, shift, index, index1, SigSize, octave, freq;
    double cos_phi, sin_phi, *value, tmp;

    if (trans == (SIGNAL *)NULL || wrd == (WRD)NULL)
	perror("GaborGetResidue(): null input!");
    if (trans[0] == (SIGNAL)NULL)
	perror("GaborGetResidue(): internal error!");

    SigSize = trans[0]->size>>1;
    cos_phi = cos((double)wrd->index->phase);
    sin_phi = sin((double)wrd->index->phase);
    shift = (int)wrd->index->position;
    octave = (int)wrd->index->octave;
    freq = (int)wrd->index->id;

    if (shift > SigSize || shift < 0)
	perror("GaborGetResidue(): illegal translation number!");
    else
	shift = SigSize - shift;

    if (octave == 0)
	{
	/*
	 * the case of dirac basis
	 */
	wrd->value = 1.0;
	trans[0]->values[(int)wrd->index->position] = 0.0;
	return;
	}
    value = trans[0]->values;
    if (octave == num_octave)
	{
/*
 * the case of Fourier basis
 */
	/*
	 * normalize the basis
	 */
	tmp = wrd->coeff*wrd->value;
	for (i=0;i<SigSize;i++)
	    {
	    index = (freq*i)%SigSize;
	    *value++ -= tmp*
		(filter[num_octave]->values[index]*cos_phi-
		filter[num_octave]->values[index+SigSize]*sin_phi);
	    }
	return;
	}
/*
 * the case of the gabor fuctions
 */
    tmp = wrd->coeff*wrd->value;
    for (i=0;i<SigSize;i++)
	{
	index = (freq*i)%SigSize;
	index1 = (shift+i)%SigSize;
	if (index1<filter[octave]->size)
	    *value++ -=  tmp*filter[octave]->values[index1]*
		(filter[num_octave]->values[index]*cos_phi-
		filter[num_octave]->values[SigSize+index]*sin_phi);
	else if (SigSize-index1<filter[octave]->size)
	    *value++ -=  tmp*filter[octave]->values[SigSize-index1]*
		(filter[num_octave]->values[index]*cos_phi-
		filter[num_octave]->values[SigSize+index]*sin_phi);
	else
	    value++;
	}
   return;
}
/*
 *
 * Get the maximum from an array provided with size, j, k, and sample rate
 * in translation
 *
 * See also complex_array_max() in complex_op.c
 *
 */
void GaborArrayMax(filter,value,size,SigSize,Log2SigSize,octave,
		freq,SampleRate_n,modula,v_r,v_i,index)
FILTER *filter;
double *value;
int size;
long SigSize;
int Log2SigSize;
int octave;
int freq;
int SampleRate_n;
double *modula;
double *v_r;
double *v_i;
int *index;
{
    int i, halfSize;
    double m, *p_r, *p_i; //, nrm=1.0;

    *modula = 0.0;
    *v_r = 0.0;
    *v_i = 0.0;
    *index = 0;
    halfSize = size>>1;
    p_r = value;
    p_i = value+size;
    for (i=0;i<size;i++)
        {
        m = (*p_r)*(*p_r)+(*p_i)*(*p_i);
        if (m>*modula)
           {
           *modula = m;
           *v_r = *p_r;
           *v_i = *p_i;
           *index = i;
           }
        p_r++;
        p_i++;
        }
}
/*
 * get the norm of the real gabor basis
 *
 * global variables used:
 *	INDEX indexG1, indexG2: indeces for inner product
 *	int cur_l, cur_h: oversample octaves for fine grid
 *	double *cur_norm:	norm of the complex gabor functions
 *	int *pnAep, *pnIep: variables for calculating the inner product
 *	double *pfCE1, *pfCE2, *pfCE3, *pfG: same as above
 *	doubel *pfC: same as above
 *	double epslonG:	the threshold for Gaussian functions
 *	int nBoundG:	the bound for Gaussina functions calculation
 *	
 * Bugs:
 *    require N is a power of 2
 */
double GaborGetGaborNorm(modula,c_r,c_i,SigSize,num_octave,
		octave,freq,translation)
double modula; /* modula square */
double c_r, c_i;
long SigSize;
int num_octave;
int octave;
long freq;
long translation;
{
    double norm, innR, innI, C;
    double cos2Phi, sin2Phi;
    double GaborGetCoeff();
    double GaborGetInnerProd();
    INDEX AllocIndex();

    if (freq==0 || freq == (SigSize>>1))
	return(1.0);

    /* set the indeces */
    if (indexG1==(INDEX)NULL)
	indexG1 = AllocIndex();
    if (indexG2==(INDEX)NULL)
	indexG2 = AllocIndex();
    indexG1->octave = indexG2->octave = (double)octave;
    indexG1->position = indexG2->position = (double)translation;
    indexG1->id = (double)freq;
    indexG2->id = (double)(SigSize-freq);
    /* calculate the norms of the complex gabor functions */
    C = GaborGetCoeff(cur_norm,octave,octave,pfC);
    /* calculate $e^{2\phi}$ */
    cos2Phi = (c_r*c_r-c_i*c_i)/modula;
    sin2Phi = 2.0*c_r*c_i/modula;

    GaborGetInnerProd(indexG1,indexG2,num_octave,cur_l,cur_h,
		0,-2*freq,nBoundG,C,
		pnAep,pnIep,
		pfCE1,pfCE2,pfCE3,pfG,
		&innR,&innI);
    norm = 2.0/(1.0+CMREAL(cos2Phi,sin2Phi,innR,innI));

    return(sqrt(norm));
}
int GaborGetSizeTime(octave,nL,nS)
int octave;
int nL;
int nS;
{
    int size;

    if (octave<nS)
	size = 1<<nL;
    else
	size = 1<<(nL-octave+nS-1);

    return(size);
}
int GaborGetSizeFreq(octave,nL,nS)
int octave;
int nL;
int nS;
{
    int size;

    if (octave>nL-nS)
	size = 1<<(nL-1);
    else
	size = 1<<(octave+nS-2);

    return(size);
}
/*
 * get the waveform
 *
 * Inputs:
 *	wrd		the information of the waveform (coeff, index, etc) (WRD)
 *	pfFilter	the filter for the decomposition (FILTER *)
 *	N		the signal size (int)
 *	L		the maximum octave (int)
 *
 * Output:
 *	sWvForm		the signal of the corresponding wave form (SIGNAL)
 */
SIGNAL GaborGetWaveForm(wrd,pfFilter,N,L)
WRD wrd;
FILTER *pfFilter;
int N;
int L;
{
    int i, octave, position, freq, index, index1;
    double cos_phi, sin_phi, *v, norm;
    SIGNAL signal, new_signal();

    octave = (int)wrd->index->octave;
    position = (int)wrd->index->position;
    freq = (int)wrd->index->id;
    norm = wrd->value;
    cos_phi = norm*cos((double)wrd->index->phase);
    sin_phi = norm*sin((double)wrd->index->phase);
    signal = new_signal(N);
    v = signal->values;
/*
 * the case of dirac
 */
    if (octave==0)
	{
	if (fabs((double)wrd->index->phase) < 1.0e-5)
	    v[position] = 1.0;
	else
	    v[position] = -1.0;
	return(signal);
	}
/*
 * the case of Fourier basis
 */
    if (octave == L)
	{
	for (i=0;i<N;i++)
	    {
	    index = (freq*i)%N;
	    *v++ = (pfFilter[L]->values[index]*cos_phi-
                pfFilter[L]->values[index+N]*sin_phi);
	    }
	return(signal);
	}
/*
 * the case of Gabor function
 */
    for (i=0;i<N;i++)
	{
	index = (freq*i)%N;
	index1 = (N-position+i)%N;
	if (index1 < pfFilter[octave]->size)
	    *v++ = pfFilter[octave]->values[index1]*
		(pfFilter[L]->values[index]*cos_phi-
		pfFilter[L]->values[index+N]*sin_phi);
	else if (N-index1<pfFilter[octave]->size)
	    *v++ = pfFilter[octave]->values[N-index1]*
		(pfFilter[L]->values[index]*cos_phi-
		pfFilter[L]->values[index+N]*sin_phi);
	else
	    *v++ = 0.0;
	}

    return(signal);
}
/*
 * Newton method to refine the maxima in the finer grid
 *
 * Inputs:
 * wrd:	selected wrd in the coarse grid (WRD)
 * trans:	transformation matrix (SIGNAL *)
 * L:		log2(N)
 * l:		subsample level in time in the fine grid (int)
 * h:		subsample level in frequency in the fine grid (int)
 * lc:		subsample level in time in the coarse grid (int)
 * hc:		subsample level in frequency in the coarse grid (int)
 *
 * Output:
 * wrd:	the selected wrd after refinement (WRD)
 */
void getMaxFrmNewton(wrd,trans,filter, MinOctave, MaxOctave, L,l,h,lc,hc)
WRD wrd;
SIGNAL *trans;
FILTER *filter;
int MaxOctave; /* Max octave used in the decomposition over gabor functions */
int MinOctave; /* Min octave used in the decomposition over gabor functions */
int L;
int l;
int h;
int lc;
int hc;
{
    long km, pm, kc, pc;
    long k, p, N;
    int jm, jc, j;
    int i, dkC, dpC, dkF, dpF, lc2, hc2, lf2, hf2;
    long badK, N2;
    double vr, vi, modula;
    double f[3][3], alpha, fm, coeff, nrm;
    double interpolation();
    double GaborGetGaborNorm();
    int getNbhd();
    void innSigWvForm();

    j = (int)wrd->index->octave;
    if (j==0 || j==L)
	return;
    k = (long)wrd->index->id;
    p = (long)wrd->index->position;
    N = 1<<L;
    N2 = N>>1;

    fm = 0.0;
    jm = j;
    km = k;
    pm = p;
    for (i=-1;i<2;i++)
	{
	jc = j+i;
	if (jc==0 || jc==L)
	    continue;
	lc2 = SSFACTOR(jc,lc);
	hc2 = SSFACTOR(L-jc,hc);
	lf2 = SSFACTOR(jc,l);
	hf2 = SSFACTOR(L-jc,h);
	dpC = STEP(jc,lc2);
	dkC = STEP(L-jc,hc2);
	dpF = STEP(jc,lf2);
	dkF = STEP(L-jc,hf2);
	/* get neighborhood */
	kc = k;
	pc = p;
	switch (getNbhd(f,trans[jc - MinOctave + 1]->values,jc,&kc,&pc,lc2,hc2,L,i))
	    {
	    case -2:
		continue;
	    case -1:
		alpha = interpolation(f,&kc,&pc,(double)dkC,(double)dpC,
			(double)dkF,(double)dpF,N,-1);
		break;
	    case 1:
		/*
		 * calculate the approximated modula
		 */
		alpha = interpolation(f,&kc,&pc,(double)dkC,(double)dpC,
			(double)dkF,(double)dpF,N,0);
		break;
	    default:
		perror("internal error! (getMaxFrmNewton())");
	    }
	badK = 1<<(L-jc);
	if ((kc>0&&kc<badK) || (kc>(N2-badK)&&kc<N2))
	    continue;
	if (alpha>fm)
	    {
	    jm = jc;
	    km = kc;
	    pm = pc;
	    fm = alpha;
	    }
	}
    if (jm==j && km==k && pm==p)
	return;
    /*
     * calculate the inner product of the R^nf and the selected wave form
     */
    innSigWvForm(trans[0],filter,jm,km,pm,L,&vr,&vi);
    modula = vr*vr+vi*vi;
    /*
     * install the selected wave form into wrd
     */
    nrm = GaborGetGaborNorm(modula,vr,vi,N,L,
			jm,km,pm);
    modula = sqrt(modula);
    coeff = modula*nrm;
    if (coeff < wrd->coeff)
	return;
    wrd->coeff = coeff;
    wrd->value = nrm;
    if (vi>0.0)
	wrd->index->phase = acos(vr/modula);
    else
	wrd->index->phase = -acos(vr/modula);
    wrd->index->octave = (double)jm;
    wrd->index->id = (double)km;
    wrd->index->position = (double)pm;

    return;
}
/*
 * get neighborhood
 */
int getNbhd(f,vector,j,k,p,l2,h2,L,flag)
double f[3][3];
double *vector;
int j;
long *k, *p;
int l2, h2;
int L;
int flag;
{
    long indxReal, indxImag;
    int Ljl2;
    long k2, p2, pw, kw, pBnd, L2, badK;
    int i, m;
    double vr, vi;
    double findMaxInNbdhd();

    L2 = 1<<(j+h2-1);
    Ljl2 = L-j+l2;
    pBnd = 1<<Ljl2;
    p2 = (*p)>>(j-l2);
    k2 = (*k)>>(L-j-h2);
    badK = 1<<h2;
    if ((k2<badK && k2>0)||(k2>(L2-badK)&&k2<L2))
	return(-2);
    if (flag==0)
	{
	indxReal = TINDXREAL(k2,p2,Ljl2);
	indxImag = TINDXIMAG(k2,p2,Ljl2);
	vr = vector[indxReal];
	vi = vector[indxImag];
	f[1][1] = vr*vr+vi*vi;
	}
    else
	{
	/* find the manxima in the nbdbood of (k2,p2) */
	f[1][1] = findMaxInNbdhd(vector,&k2,&p2,Ljl2,L2);
	if ((k2<badK && k2>0) || (k2>=(L2-badK) && k2<L2))
	    return(-2);
	}
    if (k2==badK || k2==0 || k2==(L2-badK) || k2== L2)
	{
	pw = p2-1;
	if (pw<0)
	    pw += pBnd;
	else if (pw>pBnd)
	    pw -= pBnd;
	indxReal = TINDXREAL(k2,pw,Ljl2);
	indxImag = TINDXIMAG(k2,pw,Ljl2);
	vr = vector[indxReal];
	vi = vector[indxImag];
	f[1][0] = vr*vr+vi*vi;
	pw = p2+1;
	if (pw<0)
	    pw += pBnd;
	else if (pw>pBnd)
	    pw -= pBnd;
	indxReal = TINDXREAL(k2,pw,Ljl2);
	indxImag = TINDXIMAG(k2,pw,Ljl2);
	vr = vector[indxReal];
	vi = vector[indxImag];
	f[1][2] = vr*vr+vi*vi;
	if (k2==0 || k2==L2)
	    {
	    f[1][0] /= 2.0;
	    f[1][2] /= 2.0;
	    }
	return(-1);
	}
    for (i=-1;i<2;i++)
	{
	kw = k2+i;
	for (m=-1;m<2;m++)
	    {
	    if (i==0 && m==0)
		continue;
	    pw = p2+m;
	    if (pw<0)
		pw += pBnd;
	    else if (pw>pBnd)
		pw -= pBnd;
	    indxReal = TINDXREAL(kw,pw,Ljl2);
	    indxImag = TINDXIMAG(kw,pw,Ljl2);
	    vr = vector[indxReal];
	    vi = vector[indxImag];
	    f[i+1][m+1] = vr*vr+vi*vi;
	    }
	}
    *k = k2<<(L-j-h2);
    *p = p2<<(j-l2);
    return(1);
}
/*
 * interploation
 */
double interpolation(f,k,p,dk,dp,dkF,dpF,N,flag)
double f[3][3];
long *k;
long *p;
double dk, dp; /* step sizes in the coarse grid */
double dkF, dpF; /* step sizes in the fine grid */
long N;
int flag;
{
    double fkk, fpp, fkp, fk, fp;
    double det, ktmp, ptmp;
    double fe;
    double dks, dps;

    dks = dk*dk;
    dps = dp*dp;

    if (flag<0)
	{
	fkk = 1.0;
	fkp = fk = 0.0;
	}
    else
	{
	fkk = (f[2][1]-2.0*f[1][1]+f[0][1])/dks;
	fkp = (f[2][2]-f[2][0]-f[0][2]+f[0][0])/(4.0*dk*dp);
	fk = (f[2][1]-f[0][1])/(2.0*dk);
	}
    fpp = (f[1][2]-2.0*f[1][1]+f[1][0])/dps;
    fp = (f[1][2]-f[1][0])/(2.0*dp);
    det = fkk*fpp-fkp*fkp;
/*
if (fabs(det)<1.0e-7)
    fprintf( foutput, "det=%le\n",det);
*/
    /*
     * calculate the coordinates
     */
    ktmp = (double)(*k) - (fpp*fk-fkp*fp)/det;
    ptmp = (double)(*p) - (fkk*fp-fkp*fk)/det;
    /*
     * calculate the estimated value
     */
    dk = ktmp-(double)(*k);
    dks = dk*dk;
    dp = ptmp-(double)(*p);
    dps = dp*dp;

    fe = f[1][1]+fk*dk+fp*dp+(fkk*dks+fpp*dps)/2.0+fkp*dk*dp;

    /*
     * calculate the indeces
     */
    /* map ktmp and ptmp into positive */
    while (ktmp<0)
	ktmp += (double)N;
    while (ptmp<0)
	ptmp += (double)N;
    /* map ktmp and ptmp into interger */
    *k = (long)(ktmp/dkF+.5)*(long)dkF;
    *p = (long)(ptmp/dpF+.5)*(long)dpF;
    /* map *k and *p into [0,N] */
    *k %= N;
    if (*k>(N>>1))
	*k = N-(*k);
    *p %= N;

    return(fe);
}
/*
 * inner product of a singal with a wave from`
 */
void innSigWvForm(signal,filter,j,k,p,L,vr,vi)
SIGNAL signal;
FILTER *filter;
int j;
long k;
long p;
int L;
double *vr;
double *vi;
{
    double *vs, *vf, *vCos, *vSin, tmp;
    long t, N, indxE; //, indxF;

    N = 1<<L;
    vs = signal->values;
    vf = filter[j]->values;
    vCos = filter[L]->values;
    vSin = vCos+N;
    *vr = 0.0;
    *vi = 0.0;

    for (t=0;t<=N-p-1;t++)
	{
	if (t<filter[j]->size)
	    tmp = vf[t];
	else if (N-t<filter[j]->size)
	    tmp = vf[N-t];
	else
	    tmp = 0.0;
	if (tmp==0.0)
	    continue;
	tmp *= vs[t+p];
	indxE = (k*(p+t))%N;
	*vr += tmp*vCos[indxE];
	*vi -= tmp*vSin[indxE];
	}
    for (t=1;t<=p;t++)
	{
	if (N-t<filter[j]->size)
	    tmp = vf[N-t];
	else if (t<filter[j]->size)
	    tmp = vf[t];
	else
	    tmp = 0.0;
	if (tmp==0.0)
	    continue;
	tmp *= vs[p-t];
	indxE = (k*(p-t))%N;
	*vr += tmp*vCos[indxE];
	*vi -= tmp*vSin[indxE];
	}
}
/*
 * find the maxima in the nbdhd of a given point
 */
double findMaxInNbdhd(vector,k2,p2,Ljl2,L2)
double *vector;
long *k2, *p2;
int Ljl2;
long L2;
{
   int k, p;
   long kw, pw, pBnd, pm, km;
   long indxReal, indxImag;
   double fm, f, vr, vi;

   fm = 0.0;
   km = *k2;
   pm = *p2;
   pBnd = 1<<Ljl2;
   for (k=-1;k<2;k++)
	{
	kw = (*k2)+k;
	if (kw<0 || kw>L2)
	    continue;
	for (p=-1;p<2;p++)
	    {
	    pw = ((*p2)+p+pBnd)%pBnd;
	    indxReal = TINDXREAL(kw,pw,Ljl2);
	    indxImag = TINDXIMAG(kw,pw,Ljl2);
	    vr = vector[indxReal];
	    vi = vector[indxImag];
	    f = vr*vr+vi*vi;
	    if (kw==L2)
		f /= 2.0;
	    if (f>fm)
		{
		pm = pw;
		km = kw;
		fm = f;
		}
	    }
	}
    *k2 = km;
    *p2 = pm;
    return(fm);
}
/*
 * end of ng_oper.c
 */
