#include "mpp.h"
#include "cwt1d.h"
#define BOUNDN(j,epslon,sigma) ((double)(1<<(j))*(sigma)*sqrt(-2.0*log((epslon))))
#define INDEXF(j,L,n) (abs(n)<<((L)-(j)-1))


extern double epslonG;

extern void clear_signal();

/*
 * allocation a gabor filter according to the number of octaves
 *
 * filter structure:
 *	see the manual about the gabor dictionary implementation
 *
 * Inputs:
 *	MaxOctave:	number of octaves (int)
 * Output:
 *	pfilter:	filter array (FILTER pointer)
 *
 */
FILTER *GaborAllocFilter(MaxOctave)
int MaxOctave;
{
    FILTER *pfilter;
    int j;

    if ((pfilter=(SIGNAL *)malloc((unsigned)((MaxOctave+1)*
		sizeof(struct signal))))== (SIGNAL *)NULL)
	perror("GaborAllocFilter(): mem. alloc. for pfilter failed!");

    for (j=0;j<=MaxOctave;j++)
	pfilter[j] = (FILTER)NULL;

    return(pfilter);
}
/*
 * free the gabor filter
 */
FILTER *GaborFreeFilter(pfilter,NumFilter)
FILTER *pfilter;
int  NumFilter;
{
    int j;

    if (pfilter == (FILTER *)NULL)
	return((FILTER *)NULL);
    if (NumFilter > 0)
	{
	for (j=0;j<NumFilter;j++)
	    {
	    if (pfilter[j] != (SIGNAL)NULL)
		{
		clear_signal(pfilter[j]);
		free((char *)pfilter[j]);
		}
	    }
	}
    free((char *)pfilter);
    pfilter = (FILTER *)NULL;
    return(pfilter);
}
/*
 * build a gabor filter according to the number of octaves and the signal size
 *
 * filter structre:
 * 	see the manual about the gabor dictionary implementation
 *
 * Inputs:
 *	pfilter:	filter array (FILTER pointer)
 *	MaxOctave:	number of octaves (int)
 *	SigSize:	signal size (int)
 *	sigma:		the deviation of the gaussian (double)
 *	fNorm:		the normalization factors (double *)
 * Output:
 *	pfilter:	filter array (FILTER pointer)
 *
 * Bugs:
 *	SigSize must be a power of 2
 */
FILTER *GaborBuildFilter(MaxOctave,SigSize,sigma,fNorm)
int MaxOctave;
int SigSize;
double sigma;
double **fNorm;
{
    FILTER *pfilter=(FILTER *)NULL;
    FILTER filterL;
    double boundN, boundR, boundL, B1;
    int boundI, multN, index;
    int i, j, n;
    double *value;
    double norm, factor;
    SIGNAL Gaussian2Signal(), new_signal();
    FILTER *GaborAllocFilter();
    void CreateExponential(), change_signal();
    SIGNAL new_signal(), FreeSignal();

    pfilter = GaborAllocFilter(MaxOctave);
    *fNorm = (double *)malloc(sizeof(double)*MaxOctave);
    if (*fNorm == (double *)NULL)
	perror("mem. alloc. failed!");
/*
 * calculate the non periodic Gaussian over the
 * finest grid
 */
    boundN = BOUNDN(MaxOctave-1,epslonG,sigma);
/*
    filterL = new_signal((int)boundN);
    for (n=0;n<boundN;n++)
	{
	t = (double)n/(double)(1<<(MaxOctave-1));
	filterL->values[n] = exp(-M_PI*t*t);
	}
*/
    filterL = Gaussian2Signal(0.0,boundN,
		boundN,(int)boundN+1,(double)(1<<(MaxOctave-1))*sigma);
/*
 * calculate the basic gabor functions and put it into the filter
 */
    for (j=1;j<MaxOctave;j++)
	{
	boundN = BOUNDN(j,epslonG,sigma);
	boundI = (int)(boundN/(double)SigSize+1.5);
	B1 = MIN(boundN,(double)(1<<(MaxOctave-1)));
	pfilter[j] = new_signal((int)B1+1);
	value = pfilter[j]->values;
	factor = sqrt((double)(1<<(MaxOctave-j-1)));
	for (i=0;i<=boundI;i++)
	    {
	    multN = i<<MaxOctave;
	    boundR = MIN(boundN-(double)multN,B1);
	    for (n=0;n<=(int)boundR;n++)
		{
		index = INDEXF(j,MaxOctave,n+multN);
if (index>=filterL->size)
    continue;
	        value[n] += filterL->values[index]*factor;
		}
	    if (i==0)
		continue;
	    boundR = MIN(boundN+(double)multN,B1);
	    boundL = MAX(0,(double)multN-boundN);
	    for (n=(int)boundL;n<=(int)boundR;n++)
		{
		index = INDEXF(j,MaxOctave,n-multN);
if (index>=filterL->size)
    continue;
		value[n] += filterL->values[index]*factor;
		}
	    }
	B1 = pfilter[j]->size-1;
	norm = 0.0;
	for (n=0;n<=(int)B1;n++)
	    norm += value[n]*value[n];
	norm = 2.0*norm - value[0]*value[0];
	if ((int)B1==(1<<(MaxOctave-1)))
	    norm -= value[(int)B1]*value[(int)B1];
	norm = sqrt(norm);
	(*fNorm)[j] = 1.0/norm;
	for (n=0;n<=(int)B1;n++)
	    value[n] /= norm;
	}
/*
 * store the exp(i2*pi*k/N)
 */
    CreateExponential((double)(PI2/SigSize),0.0,(double)SigSize,SigSize,
		&pfilter[MaxOctave]);

    filterL = FreeSignal(filterL);
    return(pfilter);
}
/*
 * end of gb_filter.c
 */
