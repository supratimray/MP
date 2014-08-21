/*..........................................................................*/
/*                                                                          */
/*      ------------------------------------------------------------*/
/*      (C) 1993 Copyright New York University, All Right Reserved.         */
/*	Modify by Zhifeng Zhang and Mike Orszag, 1992                       */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  signal_functions10.c                                                    */
/*                   Miscellaneous useful functions on signals:             */
/*                       Input : 1 signal                                   */
/*                       Output: 0 signal                                   */
/*                                                                          */
/****************************************************************************/
#include "mpp.h"

extern double sig_mean(); /* computes the mean of a signal */
extern double sig_variance(); /* computes the variance of a signal */

/*********************************/
/* Put a signal to zero          */
/*********************************/
sig_zero(signal)
     SIGNAL signal;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] = 0.0;
}

/*********************************/
/* Compute the absolute value of */
/* a signal.                     */
/*********************************/
sig_abs(signal)
     SIGNAL signal;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] = (double) fabs(signal->values[j]);
}
/*********************************/
/* Compute the log of a signal   */
/*********************************/
sig_log(signal)
	SIGNAL signal;
{
	int j;
	
	for (j=0; j < signal->size; j++)
	  signal->values[j] = (double) log(fabs(signal->values[j]));
}

/*
 * compute the exp of a signal
 */
sig_exp(signal)
SIGNAL signal;
{
    int j;

    if (signal==(SIGNAL)NULL)
	error("sig_exp(): null input!");
    if (signal->values==(double *)NULL)
	error("sig_exp(): null point!");

    for (j=0; j < signal->size; j++)
	signal->values[j] = (double)exp(signal->values[j]);
}
/*********************************/
/* Compute the square of a signal*/
/*********************************/
sig_square(signal)
     SIGNAL signal;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] = SQUARE(fabs(signal->values[j]));
}



/*********************************/
/* Add a number to the signal    */
/*********************************/
sig_add_num(signal, num)
     SIGNAL signal;
     double num;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] += num;
}

/*********************************/
/* Add a number to the signal    */
/*********************************/
sig_sub_num(signal, num)
     SIGNAL signal;
     double num;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] -= num;
}


/*********************************/
/* Add a number to the signal    */
/*********************************/
sig_mult_num(signal, num)
     SIGNAL signal;
     double num;
{
  int j;

  for (j = 0; j < signal->size; j++)
    signal->values[j] *= num;
}


/*********************************/
/* Add a number to the signal    */
/*********************************/
sig_div_num(signal, num)
SIGNAL signal;
double num;
{
  int j;

  if (!num)
    return;
  for (j = 0; j < signal->size; j++)
    signal->values[j] /= num;
}

/*********************************/
/* Compute the position of the   */
/* the min and the max of a      */
/* signal                        */
/*********************************/
sig_pos_min_max(signal,pmin,pmax)
     SIGNAL signal;
     double *pmin, *pmax;
{
  int j;
  double vmin, vmax;
  
  vmin = 999999;
  vmax = -999999;
  for(j=0;j<signal->size;j++)
    {
      if (vmax < signal->values[j]) {
	vmax = signal->values[j];
	*pmax = (double)j;
      }
      if (vmin > signal->values[j]) {
	vmin = signal->values[j];
	*pmin = (double)j;
      }
    }
}
/*********************************/
/* Computes the min and the max  */
/* of a signal                   */
/*********************************/
sig_min_max(signal,vmin,vmax)
     SIGNAL signal;
     double *vmin,*vmax;
{
  int j;
  
  (*vmin) = MAX_VALUE;
  (*vmax) = MIN_VALUE;
  for(j=0;j<signal->size;j++)
    {
      if ((*vmax) < signal->values[j]) (*vmax) = signal->values[j];
      if ((*vmin) > signal->values[j]) (*vmin) = signal->values[j];
    }
}
/*********************************/
/* Compute the mean of a signal  */
/*********************************/
double sig_mean(input)
     SIGNAL input;
{
  int j;
  double sum;

  sum = 0;
  for (j = 0; j < input->size; j++)
    sum += input->values[j];
  return(sum/input->size);
}
/*********************************/
/* Compute the variance of a     */
/* signal                        */
/*********************************/
double sig_variance(input)
     SIGNAL input;
{
  SIGNAL input_sq = new_struct_signal();
  double var;
  
  var = - SQUARE(sig_mean(input));
  sig_copy(input,input_sq);
  sig_square(input_sq);
  var += sig_mean(input_sq);
  delete_signal(input_sq);
  return(var);
}
/*********************************/
/* Compute the integral of a     */
/* signal                        */
/*********************************/
double integral(signal)
     SIGNAL signal;
{
  int j;
  double sum;

  sum = 0;
  for (j = 0; j < signal->size; j++)
    sum += signal->values[j] ;
  return(sum);
}
/*********************************/
/* Compute the L2 norm of a      */
/* signal                        */
/*********************************/
double sig_L2_norm_sq(input)
     SIGNAL input;
{
  SIGNAL input_sq = new_struct_signal();
  double norm, integral();
  
  sig_copy(input,input_sq);
  sig_square(input_sq);
  norm = integral(input_sq);
  delete_signal(input_sq);
  return(norm);
}
/*
 * compute the L2 square norm of a real array
 */
double farray_L2_sq_norm(value,size)
double *value;
int size;
{
   int i;
   double *v, norm=0.0;

   v = value;
   for (i=0;i<size;i++)
	{
	norm += (*v)*(*v);
	v++;
	}

   return(norm);
}
/*
 * compute the L2 norm of a real array
 */
double farray_L2_norm(value,size)
double *value;
int size;
{
   int i;
   double *v, norm=0.0;

   v = value;
   for (i=0;i<size;i++)
	{
	norm += (*v)*(*v);
	v++;
	}
   norm = (double)sqrt((double)norm);

   return(norm);
}
/*
 * find the max abs value in the signal and its position
 */
void SigAbsMx(signal,MxValue,position)
SIGNAL signal;
double *MxValue;
double *position;
{
    int i;

    if (signal == (SIGNAL)NULL)
	error("SigAbsMx(): null argument!");
    *MxValue = 0.0;
    *position = 0.0;

    for (i=0;i<signal->size;i++)
	if (fabs(signal->values[i])>fabs(*MxValue))
	    {
	    *MxValue = signal->values[i];
	    *position = (double)i;
	    }
}
