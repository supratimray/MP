/*..........................................................................*/
/*                                                                          */
/*       ------------------------Author Z. Zhang---------------------       */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  mpp.h     Basic include file for all the mpp files.       		    */
/*                   Has to be included in any new mpp file.                */
/*                                                                          */
/****************************************************************************/
#include <stdio.h>  /* Basic include files */
#include <ctype.h>
#include <string.h>
#ifndef Solaris 
	#ifndef WINDOWS
		#include <strings.h>
	#endif
#endif
#ifdef WINDOWS
	#define _USE_MATH_DEFINES
#endif

#include <math.h>
#include <malloc.h>
#ifdef sparc
#include <alloca.h>
#endif
/*
 * global constants
 */
#define NO 		0
#define YES 		1
#define TRUE		1
#define FALSE		0
//#define ERROR 		(-1)
#define STRING_SIZE 	200
#define FILTERNAME_SIZE 10
#define FILTER_SIZE	128	/* max size of a filter */
#define MAX_NUM_SIGNAL	2	/* number of signals */
#define MAX_NUM_SB	1	/* number of structure books */
#define MAX_NUM_FILTER	2	/* size of filter bank vector */
#define MAX_VALUE	1.0e38  /* max value */
#define MIN_VALUE	-1.0e38 /* min value */
#define MAX_NUM_ITERATION 100 /* number of iterations when building the book */
#define WAVELET_PATH	1
#define WHOLE_PATH	0
#define KEPT		0
#define DISCARD		1
#define OCTAVEMARK	2
#define FREQMARK	4
#define SHIFTMARK	8
#define PHASEMARK	16
#define CHIRP		32
/*
 * include file for signal, must define STRING_SIZE first
 */
#include "signals.h"    /* definition of the SIGNAL structure */
/*
 * include file for command structure, must define STRING_SIZE first
 */
//#include "commands.h"   /* definition of the COMMAND structure */
/*
 * constants for flags
 */
#define A_FLAG 		0x1
#define B_FLAG		0x2
#define C_FLAG		0x4
#define D_FLAG		0x8
#define E_FLAG		0x10
#define F_FLAG		0x20
#define G_FLAG		0x40
#define H_FLAG		0x80
#define I_FLAG		0x100
#define J_FLAG		0x200
#define K_FLAG		0x400
#define L_FLAG		0x800
#define M_FLAG		0x1000
#define N_FLAG		0x2000
#define O_FLAG		0x4000
#define P_FLAG		0x8000
#define Q_FLAG		0x10000
#define R_FLAG		0x20000
#define S_FLAG		0x40000
#define T_FLAG		0x80000
#define U_FLAG		0x100000
#define V_FLAG		0x200000
#define W_FLAG		0x400000
#define X_FLAG		0x800000
#define Y_FLAG		0x1000000
#define Z_FLAG		0x2000000
#define a_FLAG		0x4000000
#define b_FLAG		0x8000000
#define c_FLAG		0x10000000
#define d_FLAG		0x20000000
#define e_FLAG		0x40000000
#define f_FLAG		0x80000000
#define PI2		2.0*M_PI
/*
 * Useful Macros
 */
/*
 * returns YES iff inf <= x <= sup
 */
#define INRANGE(inf,x,sup)  ((inf) <= (x) && (x) <= (sup)) 
/*
 * square of a number
 */
#define SQUARE(x)  ((x)*(x))
/*
 * Max of two numbers
 */
#define MAX(x,y)  ((x) > (y) ? (x) : (y))
/*
 * Min of two numbers
 */
#define MIN(x,y)  ((x) < (y) ? (x) : (y))
/*
 * Sign of a number
 */
#define SIGN(x)  ((x) < 0 ? (-1) : (x) > 0 ? 1 : 0)
/*
 * Masking Equal operator
 */
#define MEQ( flag, mask) ( (((flag) & (mask)) == (mask)) )
/*
 * Masking Not operator
 */
#define MNOT(flag, mask) ( (((flag) & (mask)) == (mask)) ? (flag)^(mask): (flag) )
/*
 *       Filter structure
 */
typedef SIGNAL FILTER;
/*
 *       Index Structure
 */
typedef struct index{
  double id;		/* frequency of structure */
  double octave;		/* scale of structure     */
  double position;	/* position of structure  */
  double phase;		/* phase for of the structure */
} *INDEX;
/*
 *       WRD
 */
typedef struct wrd{
double coeff; /* coefficient value (real part) */
double value; /* value of the wrd */
double coeff1;
double value1;
int   status; /* either KEPT or DISCARD */
INDEX index; /* index which contains which frequency, position, and scale. */
struct wrd *next; /* next node */
} *WRD;
/*
 *       BOOK
 */
typedef struct book{
  WRD first;	   /* first element in the list */
  WRD last;	   /* last element in the list */
  int size;        /* number of structures in structure book */
  double energy;    /* energy of structure book representation */
  double sigen;     /* energy of signal */
  int id;          /* structure book volume number for I.D. purposes */
  int type;        /* structure book type */
  double smax;
  double smin;     /* maximum and minimum parameters of decomp */
  int sig_size;   /* the size of signal */
} *BOOK;

#define WAVELET 0
#define GABOR   1
#define QMF     2
#define DELTAFOURIER 3
#define NEWGABOR 4
#define FOURIER 5
#define DIRAC 6
#define WINDOWFOURIER 7

/*
 * global variables
 */

/*
 * signals
 */
extern SIGNAL signals[MAX_NUM_SIGNAL];
/*
 * structure books
 */
extern BOOK library[MAX_NUM_SB];  
/*
 * current signal
 */
extern SIGNAL cur_signal;
/*
 * current signal size
 */
extern int cur_sig_size;
/*
 * current book
 */
/* extern BOOK cur_book;  */
#define cur_book (library[Current_Book])
/*
 * previous book
 */
extern BOOK old_cur_book;
/*
 * filters
 */
/*extern SIGNAL *cur_filter, *old_cur_filter; */
extern SIGNAL *old_cur_filter;
extern int Current_Book, Old_Book;
extern int filter_type[MAX_NUM_SB];
#define cur_filter (filter[Current_Book])
#define cur_filter_type (filter_type[Current_Book])
extern SIGNAL * filter[MAX_NUM_SB];
/*
 * shift octave, subsample octave in traslation
 * and subsample octave in frequency
 */
extern int cur_shift_octave;
extern int cur_SOT;
extern int cur_SOF;
/*
 * number of filters
 */
/*extern int cur_num_filter, old_cur_num_filter; */
extern int old_cur_num_filter;
#define cur_num_filter (num_filter[Current_Book])
/*
 * number of filters
 */
extern int num_filter[MAX_NUM_SB];
/*
 * functions
 */
extern char *check_format();
extern SIGNAL *AllocFilter();
extern INDEX AllocIndex();
extern WRD AllocWord();
extern BOOK AllocBook();

extern int plot_var;

#define BW 0
#define CP 1
/*
 * input/ output
 */
extern FILE *foutput;


/*-------------------------------------------------------------------------*/
/* Function prototypes */
/*-------------------------------------------------------------------------*/

#ifdef PROTOTYPES

#include "alloc.h"
#include "int.h"
#include "book.h"
#include "sig.h"

#ifndef sun4
extern double log2(double);
#endif

#else
#ifndef sun4
extern double log2();
#endif

#endif


/*-------------------------------------------------------------------------*/

/*
 * end of mpp.h
 */
