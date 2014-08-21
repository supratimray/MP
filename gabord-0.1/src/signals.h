
/*..........................................................................*/
/*                                                                          */
/*   ------------------------Author Z. Zhang--------------------- */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/



/****************************************************************************/
/*                                                                          */
/*  signals.h        Definition of the SIGNAL structure                     */
/*                   Is included in mpp.h                            */
/*                                                                          */
/****************************************************************************/


#define SIG_SIZE 32768 /* maximum size of a signal when  (bylo 22000)
			  read from an ascii file (see file
			  signal_io.c) */



/***************************/
/*  Signal structure       */
/***************************/

typedef struct signal{
  int size_alloca;         /* size of the allocation of the 'values' field*/
  int  size;               /* size of signal */
  double shift;             /* shifting of the signal with respect to zero */
  double *values;           /* signal values */
  double scale;             /* signal scale */
  int firstp;              /* index of the first point not 
			      affected by left side effect */
  int lastp;               /* index of the last point not 
			      affected by right side effect */
  double param;             /* distance between two successive 
			      uncorrelated points */
  char name[STRING_SIZE];  /* name of the signal (not used yet) */
} *SIGNAL;


/* Functions in signal_alloc.c */

extern double *farray_malloc(); /* allocate an array of double */
extern SIGNAL new_struct_signal(); /* allocate a signal structure */
extern SIGNAL new_signal(); /* allocate a signal structure and
			       an array of double (in this signal) */







