/* ************************************************************************ */
/*		Operations on complex signals				    */
/* 									    */
/* 	The structure ofcomplex signal is SIGNAL, as defined in 	    */
/*	file signals.h . 						    */
/* 	The size of a complex signal is 2N, with N number of complex points */
/*									    */
/* 	The real part is stored in s->values[0], ..., s->values[N-1]	    */
/*	The imaginary part is  in  s->values[N], ..., s->values[2N-1]	    */
/*									    */
/* ************************************************************************ */
#include <stdio.h>
#include <math.h>
#include "mpp.h"

extern SIGNAL temporary;

/*
 *	functions defined in this file:
 *
 *	ComplexConjugate();		complex conjugate 	
 *	ComplexMultiplication();	Multiplication of 2 complex signals  
 *	ComplexMulCoeff();		Multiplication of a signal by a scalar 
 *	ComplexSubstraction();		Substraction of 2 complex signals
 *	ComplexAddition();		Addition	" 	"
 *	ComplexLinComb();		Linear Combination of 2 complex signals
 *
 *	CreateExponential();		Create function e^{-ix omega}
 *
 *	Real2Complex();   change a real signal into a complex one
 *	Complex2Real();   Extract the real part of the signal
 *
 */
/* ************************************************************************ */
/* 		Complex Conjugate of a signal				    */
/* 									    */
/*	in: SIGNAL, *SIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void ComplexConjugate(sigin, psigout)
SIGNAL sigin;
SIGNAL *psigout;
{
int index;
int HalfSize;

if ( sigin->size %2 != 0 )
	error("ComplexConjugate: size of complex signal should be even");

HalfSize = sigin->size/2;

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

/* sigout = sigin ? */
if ( *psigout != sigin )
	change_signal(*psigout, sigin->size);

for (index = 0; index <=  HalfSize-1; index ++)
	{
	(*psigout)->values[index] = sigin->values[index];
	(*psigout)->values[HalfSize + index] = -1.*sigin->values[HalfSize + index];
	}
return;
} /* end of ComplexConjugate */
/* ************************************************************************ */
/* 		Complex multiplication  of two signals			    */
/* 									    */
/*	in: SIGNAL, SIGNAL, *SIGNAL					    */
/*	out: 								    */
/* ************************************************************************ */
void ComplexMultiplication(sig1, sig2, psigout)
SIGNAL sig1, sig2;
SIGNAL *psigout;
{
int index;
int HalfSize;
double temp; /* very important temporary variable */

if ( (sig1 == NULL) || (sig2 == NULL) )
	error("ComplexMultiplication: input signal not allocated");

if ( sig1->size != sig2->size )
	error("ComplexMultiplication: size the signals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	error("ComplexMultiplication: size of complex signal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_signal(*psigout, sig1->size);

HalfSize = sig1->size/2;

for (index = 0; index <=  HalfSize-1; index ++)
	{
	/* Real part ( saved in temp ) 			*/
	/* Does not affect (*psigout)->values[index] 		*/
	/* Which is used later					*/
	/* Avoids a big problem when *psigout = sig1 or sig2 !! */

	temp = sig1->values[index] * sig2->values[index];
	temp -= sig1->values[HalfSize + index] * sig2->values[HalfSize + index];
	
	/* Imaginary part */
	(*psigout)->values[HalfSize + index] = sig1->values[index] * sig2->values[HalfSize + index] + 	sig1->values[HalfSize + index] *sig2->values[index];
	(*psigout)->values[index] = temp;
	}
return;
} /* end of ComplexMultiplication */
/* ************************************************************************ */
/* 		Complex Multiplication by a scalar of a signal		    */
/* 									    */
/*	in: SIGNAL, double, *SIGNAL					    */
/*	out: 							    */
/* ************************************************************************ */
void ComplexMulCoeff(sigin, coeff, psigout)
SIGNAL sigin;
double coeff;
SIGNAL *psigout;
{
int index;

if (sigin == NULL)
	error("ComplexMulCoeff: input signal not allocated");

if ( sigin->size %2 != 0 )
	error("ComplexMulCoeff: size of complex signal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

if(*psigout != sigin) change_signal(*psigout, sigin->size);


for (index = 0; index <=  sigin->size - 1; index ++)
	(*psigout)->values[index] = coeff * sigin->values[index];


return;
} /* end of ComplexMulCoeff */
/* ************************************************************************ */
/* 		Complex substraction  of two signals			    */
/* 									    */
/*	in: SIGNAL, SIGNAL, *SIGNAL					    */
/*	out: 								    */
/* ************************************************************************ */
SIGNAL ComplexSubstraction(sig1, sig2, psigout)
SIGNAL sig1, sig2;
SIGNAL *psigout;
{
int index;

if ( (sig1 == NULL) || (sig2 == NULL) )
	error("ComplexMultiplication: input signal not allocated");

if ( sig1->size != sig2->size )
	error("ComplexSubstraction: size the signals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	error("ComplexSubstraction: size of complex signal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_signal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = sig1->values[index] - sig2->values[index];

return;
} /* end of ComplexSubstraction */
/* ************************************************************************ */
/* 		Complex addition  of two signals			    */
/* 									    */
/*	in: SIGNAL, SIGNAL						    */
/*	out: SIGNAL							    */
/* ************************************************************************ */
SIGNAL ComplexAddition(sig1, sig2, psigout)
SIGNAL sig1, sig2;
SIGNAL *psigout;
{
int index;

if ( (sig1 == NULL) || (sig2 == NULL) )
	error("ComplexMultiplication: input signal not allocated");

if ( sig1->size != sig2->size )
	error("ComplexAddition: size the signals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	error("ComplexAddition: size of complex signal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_signal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = sig1->values[index] + sig2->values[index];

return;
} /* end of ComplexAddition */
/* ************************************************************************ */
/* 		Linear combination  of two signals			    */
/* 									    */
/*	in: double, SIGNAL, double, SIGNAL, *SIGNAL			    */
/*	out: 								    */
/* 									    */
/*	computes lambda * sig1 + mu * sig2				    */
/* ************************************************************************ */
void ComplexLinComb(lambda, sig1, mu, sig2, psigout)
double lambda;
SIGNAL sig1;
double mu;
SIGNAL sig2;
SIGNAL *psigout;
{
int index;

if ( (sig1 == NULL) || (sig2 == NULL) )
	error("ComplexMultiplication: input signal not allocated");

if ( sig1->size != sig2->size )
	error("ComplexLinComb: size the signals must be equal");

if (( sig1->size %2 != 0 ) || ( sig2->size %2 != 0 ))
	error("ComplexLinComb: size of complex signal should be even");

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

/* sigout = sig1 or sig2 ? */
if ( (*psigout != sig1) && (*psigout != sig2) )
	change_signal(*psigout, sig1->size);


for (index = 0; index <=  sig1->size-1; index ++)
	(*psigout)->values[index] = (lambda * sig1->values[index]) + (mu * sig2->values[index]);

} /* end of ComplexLinComb */
/* ************************************************************************ */
/* 			Create Exponential function			    */
/* 									    */
/*	in: double, double, double, int, *SIGNAL				    */
/*	out: 								    */
/* 									    */
/*	computes exp(i*position*omega)					    */
/* 	between min and max - scale					    */
/*	where scale = (max - min)/size					    */
/*	The length of the signal is 2*size, because it is complex	    */
/* ************************************************************************ */
void CreateExponential( position, min, max, size, psigout)
double position;
double min, max;
int size; /* size of the real signal */
SIGNAL *psigout;
{
int index;
double x, scale;

/* sigout allocated ? */
if ( *psigout == NULL ){
	*psigout = new_struct_signal(); /* allocation of sigout */
}
change_signal((*psigout), 2*size);
(*psigout)->scale = (double) (max-min)/size;
scale = (max - min)/(double)(size);


for ( index = 0; index <= size -1 ; index++ )
	{
	x = scale * (double)index + min;
	/* real part */
	(*psigout)->values[index] = (double) cos((double) (position * x) );

	/* imaginary part */
	(*psigout)->values[size + index] = (double) sin((double) (position * x) );

	}

return;
} /* end of CreateExponential */
/* ************************************************************************ */
/* 	Transformation of a real signal into a complex one		    */
/* 									    */
/*	in: SIGNAL, *SIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void Real2Complex(sigin, psigout)
SIGNAL sigin;
SIGNAL *psigout;
{
int i;

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

if ( *psigout != sigin)
	{
	change_signal(*psigout, 2 * sigin->size);

	/* real part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i] = sigin->values[i];
	/* imaginary part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i + sigin->size] = 0;
	}


if ( *psigout == sigin )
	{
	change_signal(temporary, 2 * sigin->size);

	/* real part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		temporary->values[i] = sigin->values[i];
	/* imaginary part */
	for ( i = 0; i <= sigin->size - 1; i ++ )
		temporary->values[i + sigin->size] = 0;

	sig_copy(temporary, *psigout);
	}
} /* end of Real2Complex */
/* ************************************************************************ */
/* 			Real part of a complex signal			    */
/* 									    */
/*	in: SIGNAL, *SIGNAL						    */
/*	out: 								    */
/* ************************************************************************ */
void Complex2Real(sigin, psigout)
SIGNAL sigin;
SIGNAL *psigout;
{
int i;
int Size;

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

Size = sigin->size / 2;

if ( *psigout != sigin)
	{
		change_signal(*psigout, Size);

	/* real part */
	for ( i = 0; i <= Size - 1; i ++ )
		(*psigout)->values[i] = sigin->values[i];
	}

if ( *psigout == sigin )
	{
	change_signal(temporary,  Size);

	/* real part */
	for ( i = 0; i <= Size - 1; i ++ )
		temporary->values[i] = sigin->values[i];

	sig_copy(temporary, *psigout);
	}
} /* end of Complex2Real */
/*
 * take the imaginary part of the complex signal
 */
void Complex2Imaginary(sigin, psigout)
SIGNAL sigin;
SIGNAL *psigout;
{
int i;
int Size;

/* sigout allocated ? */
if ( *psigout == NULL )
	*psigout = new_struct_signal(); /* allocation of sigout */

Size = sigin->size / 2;

if ( *psigout != sigin)
	{
		change_signal(*psigout, Size);

	/* real part */
	for ( i = Size; i <= sigin->size - 1; i ++ )
		(*psigout)->values[i-Size] = sigin->values[i];
	}

if ( *psigout == sigin )
	{
	change_signal(temporary,  Size);

	/* real part */
	for ( i = Size; i <= sigin->size - 1; i ++ )
		temporary->values[i-Size] = sigin->values[i];

	sig_copy(temporary, *psigout);
	}
} /* end of Complex2Real */
/*
 * taking complex modulation from a complex
 * signal
 */
void ComplexModulation(sigin,sigout)
SIGNAL sigin;
SIGNAL sigout;
{
    int i, size;
    if (sigin == (SIGNAL)NULL || sigout == (SIGNAL)NULL)
	error("ComplexModulation(): null argument!");

    size = sigin->size/2;

    if (sigout->size != size)
	change_signal(sigout,size);

    for (i=0;i<size;i++)
	sigout->values[i] = sqrt((double)(sigin->values[i]*sigin->values[i]+
		sigin->values[i+size]*sigin->values[i+size]));
}
/*
 * multiply a complex signal by a complex number
 */
void ComplexMulCN(signal,v_real,v_imag)
SIGNAL signal;
double v_real;
double v_imag;
{
    int i, size;
    double temp;

    if (signal == (SIGNAL)NULL)
	error("ComplexMulCN(): null input!");

    size = signal->size/2;
    for (i=0;i<size;i++)
	{
	temp = v_real*signal->values[i] -
			v_imag*signal->values[i+size];
	signal->values[i+size] = v_imag*signal->values[i] +
			v_real*signal->values[i+size];
	signal->values[i] = temp;
	}
}
/*
 * searching the maximum modula from a complex array
 *
 * The first half of the array contains the real part and the second
 * half of the array contains the imaginary part
 *
 * Inputs:
 *	value		the complex array (double *)
 *	size		size of the complex array (int)
 *	size_op		size for the operation (int)
 *
 * Outputs:
 *	modula		the maximum modula square found (double *)
 *	v_r		the real part of the maximum modula found (double *)
 *	v_i		the imaginary part of the maximum modula 
 *			found (double *)
 *	index		the index of the maximum modula (int)
 *
 */
void complex_array_max(value,size,size_op,modula,v_r,v_i,index)
double *value;
int size;
int size_op;
double *modula;
double *v_r;
double *v_i;
int *index;
{
    int i;
    double m, *p_r, *p_i;

    *modula = 0.0;
    *v_r = 0.0;
    *v_i = 0.0;
    *index = 0;
    p_r = value;
    p_i = value+size;
    for (i=0;i<size_op;i++)
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
 * end of complex_op.c
 */
