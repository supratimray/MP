
#include "mpp.h"
#include "cwt1d.h"
/*
 * decomposition for the gabor library
 *
 * Inputs:
 *	trans			Stored the result of the transformation
 *				(SIGNAL *)
 *	signal			Input (complex) signal to be decomposed (SIGNAL)
 *	filter			Stored the basic gabor functions (FILTER *)
 *	SubsampleOctaveTime	The octave number that begin to subsample
 *				for the translateion (int)
 *	SubsampleOctaveFreq	The octave number that begin to subsample
 *				for the frequency (int) (not used yet)
 *	MaxOctave		Number of octaves (int)
 *	ShiftOctave		The octave number that begin to use
 *				the seconf formula to calculate the
 *				transformation (int)
 * Outputs:
 *	trans			Stored the result of the transformation
 *				(SIGNAL *)
 * Bugs:
 *	1. This function handle a complex signal;
 *	2. This function does not check the
 *	   validity of 'SubsampleOctaveTime','SubsampleOctaveFreq', and 
 *	   'ShifOctave'.(for example, they must all less than MaxOctave
 *	   and consistant with the filter;
 *	3. Require 'signal->size' is a power of 2.
 */
SIGNAL *GaborDecomp(trans,signal, filter, SubsampleOctaveTime,
		SubsampleOctaveFreq, MinOctave, MaxOctave, ShiftOctave)
SIGNAL *trans;
SIGNAL signal;
FILTER *filter;
int SubsampleOctaveTime;
int SubsampleOctaveFreq;
int MinOctave, MaxOctave;
int ShiftOctave;
{
    int i, j, k, k_length, SigSize, DoubleSigSize, shift, size, SampleRate_t=1;
    int HalfSigSize, index;
    int LnSigSize;
    int n, n_length, n_2_length, SampleRate_f;
    SIGNAL sigtmp=(SIGNAL)NULL;
    double *value_r, *value_i, *value1, *value2, *value3, *value4, *value5;
    double *sig_value_r, *sig_value_i, SqrtSigSize;
    double *darray_malloc();
    SIGNAL new_signal(), *GaborAllocFilter();
    int IFFT(), FFT();  /* fft.c */
    int delete_signal(); /* sig_alloc.c */
    void change_signal(), sig_range_copy(), farray_translate();
    void GaborSubsample();

    if (trans == (SIGNAL *)NULL)
	trans = GaborAllocFilter(MaxOctave- MinOctave + 3);
    if (signal == (SIGNAL)NULL || signal->size == 0)
	perror("Empty signal!");

    DoubleSigSize = signal->size;
    SigSize = DoubleSigSize>>1;
    HalfSigSize = SigSize>>1;
    SqrtSigSize = (double)sqrt((double)SigSize);
    LnSigSize = (int)log2((double)SigSize);
/*
 * put the original signal into trans[0]
 */
    trans[0] = signal;
/*
 * put the Fourier transform of the original signal into trans[MaxOctave - MinOctave + 2]
 */
    if (trans[MaxOctave - MinOctave + 2] == (SIGNAL)NULL)
	trans[MaxOctave - MinOctave + 2] = new_signal(DoubleSigSize);
    else if (trans[MaxOctave - MinOctave + 2]->size != DoubleSigSize)
	change_signal(trans[MaxOctave - MinOctave + 1],DoubleSigSize);
    value1 = trans[MaxOctave - MinOctave + 2]->values;
    value2 = trans[MaxOctave - MinOctave + 2]->values+SigSize;
    value3 = trans[0]->values;
    value4 = trans[0]->values + SigSize;
    for (i=0;i<SigSize;i++)
	{
	*value1++ = *value3++;
	*value2++ = *value4++;
	}
    FFT(trans[MaxOctave - MinOctave + 2],&trans[MaxOctave - MinOctave + 2]);
/*
 * calculate <f,g_j,k,n> using the first formula.
 */
    if (sigtmp == (SIGNAL)NULL)
	sigtmp = new_signal(DoubleSigSize);
    sig_value_r = darray_malloc(SigSize);
    sig_value_i = darray_malloc(SigSize);
    for (j= MinOctave;j<ShiftOctave;j++)
	{
	/*
	 * calculate the length for k
	 */
	if (j>(LnSigSize-SubsampleOctaveFreq))
	    {
	    k_length = 1<<(LnSigSize-1);
	    SampleRate_f = 1;
	    }
	else
	    {
	    k_length = 1<<(SubsampleOctaveFreq+j-2);
	    SampleRate_f = 1<<(LnSigSize-SubsampleOctaveFreq-j+1);
	    }
	/*
	 * calculate the length for n
	 */
	if (j<SubsampleOctaveTime)
	    {
	    SampleRate_t = 1;
	    n_length = SigSize;
	    }
	else
	    {
	    SampleRate_t = 1<<(j-SubsampleOctaveTime+1);
	    n_length = SigSize>>(j-SubsampleOctaveTime+1);
	    }
	n_2_length = n_length<<1;
	/*
	 * calculate the size for the array <f,g_j,k,n> for j fixed
	 */
	//	size = k_length*n_2_length+n_2_length;
	size = k_length*n_2_length+n_2_length+2;  //  PJF 

	if (trans[j - MinOctave + 1] == (SIGNAL)NULL)
	    trans[j - MinOctave + 1] = new_signal(size);
	else if (trans[j - MinOctave + 1]->size != size)
	    change_signal(trans[j - MinOctave + 1],size);
	value_r = trans[j - MinOctave + 1]->values;
	for (k=0;k<=k_length;k++)
	    {
	    value1 = trans[MaxOctave - MinOctave + 2]->values;
	    value2 = trans[MaxOctave - MinOctave + 2]->values+SigSize;
	    /*
	     * calculate f^(w+2^(L-L1-j+1)*k)
	     */
	    shift = -SampleRate_f*k;
	    farray_translate(value1,sig_value_r,SigSize,shift);
	    farray_translate(value2,sig_value_i,SigSize,shift);
	    /*
	     * calculate ^g_j(p)*^f(p-2^(L-j)k)
	     */
	    value1 = sigtmp->values;
	    value2 = sigtmp->values+SigSize;
	    value3 = filter[LnSigSize-j]->values;
	    /* the Fourier transform of g_j/sqrt(SigSize) */
	    value4 = sig_value_r;
	    value5 = sig_value_i;
	    value1[0] = value3[0]*value4[0]*SqrtSigSize;
	    value2[0] = value3[0]*value5[0]*SqrtSigSize;
	    for (i=1;i<HalfSigSize;i++)
		{
		if (i<filter[LnSigSize-j]->size)
		    {
		    value1[i] = value3[i]*value4[i]*SqrtSigSize;
		    value2[i] = value3[i]*value5[i]*SqrtSigSize;
		    value1[SigSize-i] = value3[i]*value4[SigSize-i]*
						SqrtSigSize;
		    value2[SigSize-i] = value3[i]*value5[SigSize-i]*
						SqrtSigSize;
		    }
		else
		    {
		    value1[i] = value2[i] = 0.0;
		    value1[SigSize-i] = value2[SigSize-i] = 0.0;
		    }
		} /* end of i loop */
	    if (filter[LnSigSize-j]->size>HalfSigSize)
		{
		value1[HalfSigSize] = value3[HalfSigSize]*value4[HalfSigSize]*
					SqrtSigSize;
		value2[HalfSigSize] = value3[HalfSigSize]*value5[HalfSigSize]*
					SqrtSigSize;
		}
	    else
		value1[HalfSigSize] = value2[HalfSigSize] = 0.0;
	    if (j>=SubsampleOctaveTime)
		/*
		 * calculate the corresponding Fourier transform
		 * for the subsample signal {<f,g_j,k,n>} in n with
		 * rate of 2^(j-l)
		 */
		{
		for (i= SubsampleOctaveTime;i<=j;i++)
		    GaborSubsample(sigtmp);
		}
	    IFFT(sigtmp,&sigtmp);
	    value1 = sigtmp->values;
	    value2 = sigtmp->values+n_length;
	    value_i = value_r+n_length;
	    for (i=0;i<n_length;i++)
		{
/*
		*value_r++ = sqrtsr*(*value1++)/(double)SigSize;
		*value_i++ = sqrtsr*(*value2++)/(double)SigSize;
*/
		*value_r++ = *value1++;
		*value_i++ = *value2++;
		} /* end of i loop */
/*
 * reset the size of the sigtmp (changed when doing the subsampling
 */
	    sigtmp->size = DoubleSigSize;
	    value_r = value_i;
	    } /* end of k loop */
	} /* end of loop j */
/*
 * calculate <f,g_j,k,n> using the second formula
 */
	/*
	 * 2N*2^(l-2)
	 */
    for (j=ShiftOctave;j<=MaxOctave;j++)
	{
	/*
	 * calculate the length for k and subsample rate for k
	 */
	if (j>(LnSigSize-SubsampleOctaveFreq))
	    {
	    k_length = (1<<LnSigSize)-1;
	    SampleRate_f = 1;
	    }
	else
	    {
	    k_length = 1<<(SubsampleOctaveFreq+j-2);
	    SampleRate_f = 1<<(LnSigSize-SubsampleOctaveFreq-j+1);
	    }
	/*
	 * calculate the length and subsample rate for n
	 */
	if (j<SubsampleOctaveTime)
	    {
	    SampleRate_t = 1;
	    n_length = SigSize;
	    }
	else
	    {
	    SampleRate_t = 1<<(j-SubsampleOctaveTime+1);
	    n_length = SigSize>>(j-SubsampleOctaveTime+1);
	    }
	n_2_length = n_length<<1;
	/*
	 * calculate the size for the array <f,g_j,k,n> for j fixed
	 */
	//	size = k_length*n_2_length+n_2_length;
	size = k_length*n_2_length+n_2_length+2;  // PJF

	if (trans[j - MinOctave + 1] == (SIGNAL)NULL)
	    trans[j - MinOctave + 1] = new_signal(size);
	else if (trans[j - MinOctave + 1]->size != size)
	    change_signal(trans[j - MinOctave + 1],size);
	for (n=0;n<n_length;n++)
	    {
	    /*
	     * calculate g_j(p-n)
	     */
	    shift = n*SampleRate_t;
/*
	    value1 = filter[j]->values;
	    farray_translate(value1,sig_value_r,SigSize,shift);
*/
	    /*
	     * calculate f(p)g_j(p-n)
	     *
	     * sigtmp->values[i] = trans[0]->values[i]*sig_value_r[i]
	     * sigtmp->values[i+SigSize] = trans[0]->values[i+SigSize]*
	     *					sig_value_r[i];
	     */
	    value1 = sigtmp->values;
	    value2 = value1 + SigSize;
	    value3 = trans[0]->values;
	    value4 = value3 + SigSize;
	    value5 = filter[j]->values;
	    for (i=0;i<SigSize;i++)
		{
		index = (SigSize+i-shift)%SigSize;
		if (index<filter[j]->size && index != 0 
			&& index != HalfSigSize)
		    {
		    value1[i] = value3[i]*value5[index];
		    value2[i] = value4[i]*value5[index];
		    }
		else if (SigSize-index < filter[j]->size && index != 0
				&& index != HalfSigSize)
		    {
		    value1[i] = value3[i]*value5[SigSize-index];
		    value2[i] = value4[i]*value5[SigSize-index];
		    }
		else if (index==0)
		    {
		    value1[i] = value3[i]*value5[0];
		    value2[i] = value4[i]*value5[0];
		    }
		else if (index==HalfSigSize)
		    {
		    if (filter[j]->size>HalfSigSize)
			{
			value1[i] = value3[i]*value5[HalfSigSize];
			value2[i] = value4[i]*value5[HalfSigSize];
			}
		    else
			value1[i] = value2[i] = 0.0;
		    }
		else
		    value1[i] = value2[i] = 0.0;
		}
	    /*
	     * calculate the ^(f*g_j,0,n)(2^(L-L1-j+1)*k)
	     */
	    /*
	     * calculate the corresponding f*g_j,0,n for subsampling
	     * in the frequency domain
	     */
	    for (i=j;i<=LnSigSize-SubsampleOctaveFreq;i++)
		GaborSubsample(sigtmp);
	    /*
	     * Fourier transform of f*g_j,0,n
	     */
	    FFT(sigtmp,&sigtmp);
	    /*
	     *
	     * put the results into trans[j - MinOctave + 
	     * 1]->values[k*n_length+n]
	     *    trans[j]->values[k*n_2_length+n] = sigtmp->values[k]*
	     *				SampleRate_f;
	     *    trans[j - MinOctave + 1]->values[k*n_2_length+n_length+n] = 
	     *			sigtmp->values[k+n_length]*SampleRate_f;
	     *
	     * SampleRate_f is the normalization factor
	     *
	     */
	    value1 = sigtmp->values;
	    value2 = value1 + (sigtmp->size>>1);
	    value_r = trans[j - MinOctave + 1]->values + n;
	    for (k=0;k<=k_length;k++)
		{
		value_i = value_r + n_length;
		*value_r = (*value1++)*(double)SampleRate_f;
		value_r += n_2_length;
		*value_i = (*value2++)*(double)SampleRate_f;
		}
	    sigtmp->size = DoubleSigSize;
	    } /* end of loop n */
	} /* end of loop j */
    /*
     * normalize the Fourier transform of the orignal signal
     */
    value1 = trans[MaxOctave - MinOctave + 2]->values;
    value2 = value1 + SigSize;
    for (i=0;i<SigSize;i++)
	{
	*value1++ /= SqrtSigSize;
	*value2++ /= SqrtSigSize;
	}
/*
 * free the local memories
 */
    free((char *)sig_value_r);
    free((char *)sig_value_i);
    delete_signal(sigtmp);

    return(trans);
}/* end of GaborDecomp() */
/*
 * calculate the Fourier transform a subsample signal at rate 2
 *
 * Inputs:
 *	signal		the Fourier transform of the original signal
 *
 * Output:
 *	signal		the Fourier transform of the subsampled signal
 *			at rate 2
 *
 * Bugs:
 *	require the signal size is a power of 2
 */
void GaborSubsample(signal)
SIGNAL signal;
{
    int newsize, new_sig_size, i;
    double *value1, *value2, *value3, *value4;

    if (signal == (SIGNAL)NULL)
	perror("GaborSubsample(): null input!");
    if (signal->size <= 0)
	return;

    newsize = signal->size>>1;
    new_sig_size = newsize>>1;
/*
 * change the signal size, but dos not change the
 * allocation size
 */
    signal->size = newsize;
/*
 * get the corresponding pointers for faster computation
 */
    value1 = signal->values;
    value2 = signal->values+new_sig_size;
    value3 = signal->values+newsize;
    value4 = value3+new_sig_size;
/*
 * calculate
 *     the real part
 *	f(p)=((f(p)+f(p+N/2))/2.0
 *     the imaginary part
 *      f(P+N/2) = ((f(p+N)+f(p+3N/2))/2.0
 */ 
    for (i=0;i<new_sig_size;i++)
	{
	*value1 += *value2;
	*value1++ /= 2.0;
	*value2++ = ((*value3++)+(*value4++))/2.0;
	}
}
/*
 * end of ng_decomp.c
 */
