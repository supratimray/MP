#include "mpp.h"
#include "cwt1d.h"
#include <sys/types.h>
#ifdef WINDOWS
	#include <WinSock2.h>
#else
	#include <sys/times.h>
	#include <sys/param.h>
#endif
/* real part of complex multiplication */
#define CMREAL(a1,b1,a2,b2)	((a1)*(a2)-(b1)*(b2))
/* imaginary part of complex multiplication */
#define CMIMAGINARY(a1,b1,a2,b2)	((a1)*(b2)+(a2)*(b1))
/*
 *
 * Build book from a gabor transform
 *
 * Inputs:
 *	trans		gabor transform (SIGNAL *)
 *	filter		the basic gabor fuctions, used for updating
 *			(FILTER *)
 * 	MinOctave	Minimum of the octave decomposition on gabor functions
 *	MaxOctave	Maximum of the octave decomposition on gabor functions
 *	threshold	precision that stop the loop (double)
 *	LnSize		Log of the signal size
 *	sig_norm	the L2 norm of the original signal (double)
 *	num_iter	current number of iterations (int)
 *	max_num_iter	number of iterations allowed in the loop
 *	ShiftOctave	the octave that begins to use the second formula
 *			compute the projection
 *	SubsampleOctaveTime	the octave that begins to subsample in
 *				translation
 *	SubsampleOctaveFreq	the octave that begins to subsample in
 *				frequency
 *
 * Outputs:
 *	book		book that stores the results (BOOK)
 *
 */
void GaborBuildBook(trans,filter,book, MinOctave, MaxOctave, LnSize,
		SubsampleOctaveTime,SubsampleOctaveFreq,
		l,h,pnIep,pfC,pfB,pfG,pfCE1,
		pfFilterNorm)
SIGNAL 	*trans;
FILTER	*filter;
BOOK	book;
int	MinOctave, MaxOctave; /* Min and Max octave for the decomposition on Gabor functions */
int	LnSize;
int	SubsampleOctaveTime;
int	SubsampleOctaveFreq;
int	l, h;	/* oversampling octaves for the fine grid */
int	*pnIep;
double	*pfC, *pfB;
double	*pfG;
double	*pfCE1;
double	*pfFilterNorm;
{

    extern void GaborGetResidue();
    WRD wrd;

    WRD GaborGetMaxFrmTrans();
    SIGNAL *GaborDecomp();
    extern void BookAppend(), UpdateGabor(), GaborUpdateFourier();
    void getMaxFrmNewton();

    if (trans == (SIGNAL *)NULL || book == (BOOK)NULL)
	perror("GaborBuildBook(): null arguments!");
	
/* get the maximum from the trans and put it into wrd*/

    wrd = GaborGetMaxFrmTrans(trans,filter, MinOctave, MaxOctave, LnSize,SubsampleOctaveTime,
		SubsampleOctaveFreq);

    if (wrd->index->octave != 0.0 
	&& wrd->index->octave != (double)LnSize
	&& (l>SubsampleOctaveTime || h>SubsampleOctaveFreq))
	getMaxFrmNewton(wrd,trans,filter, MinOctave, MaxOctave, LnSize,l-1,h-1,
		SubsampleOctaveTime-1,
		SubsampleOctaveFreq-1);

/* update the transform */
 
    GaborGetResidue(trans,filter,wrd,LnSize);

    UpdateGabor(trans,wrd, MinOctave, MaxOctave, LnSize,
		SubsampleOctaveTime-1,SubsampleOctaveFreq-1,l-1,h-1,
		pfG,pfC,pfFilterNorm,filter,pfB,pfCE1,pnIep);

   GaborUpdateFourier(trans,filter,wrd, MinOctave, MaxOctave, LnSize);		

/* put the wrd into book */

    BookAppend(book,wrd);
    return;
}


////// update the Fourier basis
 
void GaborUpdateFourier(trans,filter,wrd, MinOctave, MaxOctave, L)
SIGNAL *trans;
FILTER *filter;
WRD wrd;
int	MinOctave, MaxOctave; /* Min and Max octave for the decomposition on Gabor functions */
int L;
{
    double *v1r, *v1i, *v2r, *v2i, *vg;
    double coeff, tmp, phi1, cosPhi, sinPhi;
    long j1, k1, p1;
    long index1, index2;
    long k2, N;

    N = (trans[0]->size)>>1;
    j1 = wrd->index->octave;
    k1 = wrd->index->id;
    p1 = wrd->index->position;
    phi1 = wrd->index->phase;
    coeff = wrd->coeff*wrd->value;
    v1r = trans[MaxOctave - MinOctave + 2]->values;
    v1i = v1r+N;
    if (j1==0)
	{
	/* the case of dirac */
	v2r = filter[L]->values;
	v2i = v2r+N;
	if (fabs(phi1)>M_PI_2)
	    coeff = -coeff/sqrt((double)N);
	else
	    coeff /= sqrt((double)N);
	for (k2=0;k2<N;k2++)
	    {
	    index1 = (k2*p1)%N;
	    *v1r++ -= coeff*v2r[index1];
	    *v1i++ += coeff*v2i[index1];
	    }
	}
    else if (j1==L)
	{
	/* the case of Fourier */
	v1r[k1] = 0.0;
	v1i[k1] = 0.0;
	if (k1!=0)
	    {
	    v1r[N-k1] = 0.0;
	    v1i[N-k1] = 0.0;
	    }
	}
    else
	{
	/* the case of gabor */
	cosPhi = cos(phi1);
	sinPhi = sin(phi1);
	coeff /= 2.0;
	v2r = filter[L]->values;
	v2i = v2r+N;
	vg = filter[L-j1]->values;
	for (k2=0;k2<N;k2++)
	    {
	    index1 = k2-k1+N;
	    index2 = (index1*p1)%N;
	    index1 %= N;
	    if (index1<filter[L-j1]->size)
		tmp = vg[index1];
	    else if ((N-index1)<filter[L-j1]->size)
		tmp = vg[N-index1];
	    else
		tmp = 0.0;
	    if (tmp!=0.0)
		{
		tmp *= coeff;
		*v1r -= 
		tmp*CMREAL(v2r[index2],-v2i[index2],cosPhi,sinPhi);
		*v1i -= 
		tmp*CMIMAGINARY(v2r[index2],-v2i[index2],cosPhi,sinPhi);
		}
	    index1 = k2+k1;
	    index2 = (index1*p1)%N;
	    index1 %= N;
	    if (index1<filter[L-j1]->size)
		tmp = vg[index1];
	    else if ((N-index1)<filter[L-j1]->size)
		tmp = vg[N-index1];
	    else
		tmp = 0.0;
	    if (tmp!=0.0)
		{
		tmp *= coeff;
		*v1r -= 
		tmp*CMREAL(v2r[index2],-v2i[index2],cosPhi,-sinPhi);
		*v1i -= 
		tmp*CMIMAGINARY(v2r[index2],-v2i[index2],cosPhi,-sinPhi);
		}
	    v1r++;
	    v1i++;
	    }
	}
}
/*
 * end of ng_buildbook.c
 */
