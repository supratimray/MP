/*******************************************************************/
/*  MPP signal processing program. 			   */
/* (C) 1993 Copyright New York University, All Right Reserved.			   */
/*							           */
/*  ------------------------------------------------------------  */
/*  Francois Bergeaud, Mike Orszag.                                */
/*******************************************************************/


/****************************************************************************/
/*                                                                          */
/*  signal_alloc.c   Functions which deal with the dynamical                */
/*                   allocation of memory for SIGNAL's                      */
/*                                                                          */
/****************************************************************************/



#include "mpp.h"


/************************************/
/* allocation of an array of doubles */
/************************************/
double *darray_malloc(size)
int size;
{
  double *ptr;
  //  if(ptr = (double *) malloc(sizeof(double)*size)) return(ptr);
  if(ptr = (double *) calloc(size,sizeof(double))) return(ptr);    // PJF
  else {
    fprintf(stderr,"Can't find enough memory for %d doubles\n",size);
    error("");
  }
}
  
/* there is also a macro using "alloca" instead of "malloc" defined
   in mpp.h and called farray_alloca but it does NOT test wether
   the allocation succeeded or not */




/*************************************************/
/* Create a new signal structure and returns it  */
/*************************************************/
SIGNAL new_struct_signal()
{
  SIGNAL signal;
  void init_signal();

  if(!(signal = (SIGNAL) (malloc(sizeof(struct signal)))))
    error("Mem. alloc for SIGNAL failed\n");

  init_signal(signal);

  return (signal);
}

void init_signal(signal)
SIGNAL signal;
{
  if (signal == (SIGNAL)NULL)
	error("init_signal(): null argument!");

  signal->values = (double *)NULL;
  signal->size_alloca = 0;
  signal->size = 0;
  signal->name[0] = '\0';
  signal->scale = 1.;
  signal->shift = 0.;
  signal->firstp = 0;
  signal->lastp = 0;
  signal->param = 1.;
}


/*************************************************/
/* Desallocate the whole 'signal' structure      */
/*************************************************/
delete_signal(signal)
SIGNAL signal;
{
  if (signal)
    {
    if (signal->values) free((char *)signal->values);
	free((char *)signal);
    }
}


/*************************************************/
/* Create a signal structure and an array of     */
/* double of size 'size'. The array is put in the */
/* 'values' field of the signal.                 */
/* It returns the signal                         */
/*************************************************/
SIGNAL new_signal(size)
int size;
{
  int i;

  SIGNAL signal = new_struct_signal();
  signal->values = darray_malloc(size);
  signal->size = size;
  signal->size_alloca = size;
  signal->lastp = size-1;
  signal->param = 1.;
  for (i=0;i<size;i++)
	signal->values[i] = 0.0;
  return(signal);
}


/*************************************************/
/* Initialization of 'signal' and desallocation  */
/* of the array of double signal->values          */
/*************************************************/
void clear_signal(signal)
     SIGNAL signal;
{
   void init_signal();

/*   if (signal == (SIGNAL)NULL)
	error("clear_signal(): null argument!"); */

   if (signal->values)
	{
	free((char *)signal->values);
	init_signal(signal);
	}
}
    

/*************************************************/
/* Return YES if signal->values is not NULL      */
/*************************************************/
int is_empty(signal)
     SIGNAL signal;
{
  if (signal->values) {return(NO);} else return(YES);
}


/*************************************************/
/* Very important procedure which has to be      */
/* called each time one needs to store doubles    */
/* in a signal.                                  */
/* If 'signal' has already a size > 'size'       */
/* it doesn't do any memory allocation.          */
/* If not then it first desallocates             */
/* signal->values and replaces it by an array    */
/* of size 'size.                                */
/* In all cases, it initializes the fields       */
/* of 'signal'.                                  */
/*************************************************/
void change_signal(signal,size)
     SIGNAL signal;
     int size;
{
  int i;

  if (signal->size_alloca < size)
    {
      clear_signal(signal);
      signal->values = darray_malloc(size);
      signal->size_alloca = size;
      if (signal->values == NULL)
	printf ("malloc(%d) failed!\n",size);
    }

  signal->size = size;
  signal->scale = 1.;
  signal->shift = 0.;
  signal->firstp = 0;
  signal->lastp = size-1;
  signal->param = 1.;
  for(i = 0; i < size; i++)
    signal->values[i] = 0.0;
}
/*
 * free SIGNAL structure and return NULL
 */
SIGNAL FreeSignal(signal)
SIGNAL signal;
{
    delete_signal(signal);
    return((SIGNAL)NULL);
}
