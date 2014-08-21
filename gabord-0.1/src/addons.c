#include "mpp.h"
//extern BOOK library[MAX_NUM_SB];

extern void change_signal();

/*
 * calculate a signal created by the L2 normalized gaussian
 *
 * Inputs:
 * 	min_t	minimum value for t (double)
 *	max_t	maximum value for t (double)
 *	size_g	number of points in t need to be calculated (int)
 *	size_s  size of the signal need to be created (int),
 *	sigma	the deviation for the gaussian
 * 
 * Remark:
 *	if size_s > size_g, then the signal filled with zero for
 *	those points that index number greater than size_g
 *
 * Bugs:
 *	size_s must be greater than size_g, if not return null pointer
 *
 */
SIGNAL Gaussian2Signal(min_t,max_t,size_g,size_s,sigma)
double min_t;
double max_t;
double size_g;
int size_s;
double sigma;
{
    SIGNAL signal=(SIGNAL)NULL;
    double t, scale;
    int i;
    double gaussianL2();
    SIGNAL new_signal();

    if (size_s < size_g)
	return((SIGNAL)NULL);

    signal = new_signal(size_s);
    scale = (double)(max_t-min_t)/(double)size_g;
    signal->scale = (double)scale;

    t = (double)min_t;

    for (i=0;i<=(int)size_g;i++)
	{
	signal->values[i] = gaussianL2(t,(double)sigma);
	t += scale;
	}

    return(signal);
}


/**************************************/
/* Copy a signal in another           */
/**************************************/

/* copy an array of double ('input') of 
   size 'sigsize' in another ('output') */
/*--------------------------------------------------------------------------*/
farray_copy(input, sigsize, output)
     double *input;
     int sigsize;
     double *output;
{
  double *to, *from;

  for (to=output,from = input; from < input+sigsize; to++, from++)
		*to = *from;
}

/*--------------------------------------------------------------------------*/
/*
 * tranlate an array of double
 */
/*--------------------------------------------------------------------------*/
void farray_translate(f_in,f_out,size,shift)
double *f_in;
double *f_out;
int size;
int shift;
{
    int i, i0, i1;
    double *value;

    if (f_in == (double *)NULL || f_out == (double *)NULL)
	perror("farray_translate(): null input!");

    while (shift < 0)
	shift += size;
    while (shift-size>0)
	shift -= size;
    i0 = size-shift;
    i1 = i0+size;
    value = f_out;
    for (i=i0;i<i1;i++)
	*value++ = f_in[i%size];
}

/*--------------------------------------------------------------------------*/
/* Signal Copy */
/*--------------------------------------------------------------------------*/
sig_copy(input,output)
     SIGNAL input,output;
{
  int i;

  if (input == (SIGNAL)NULL || output == (SIGNAL)NULL)
	perror("sig_copy(): null input!");

  change_signal(output,input->size);

  for(i=0;i<input->size;i++)
    output->values[i] = input->values[i];
  output->scale = input->scale;
  output->shift = input->shift;
  output->firstp = input->firstp;
  output->lastp = input->lastp;
  output->param = input->param;
}

/*--------------------------------------------------------------------------*/
/*
 * signal copy from range min to max
 */
/*--------------------------------------------------------------------------*/
sig_range_copy(input,output,min,max)
SIGNAL input;
SIGNAL output;
int min;
int max;
{
  int i, size;
  double *value;

  if (input == (SIGNAL)NULL || output == (SIGNAL)NULL)
	perror("sig_range_copy(): null input!");
  if (min<0 || max<min || max>input->size)
	perror("sig_range_copy(): illegal argument!");

  size = max-min;
  if (output->size != size)
	change_signal(output,size);

  value = output->values;
  for (i=min;i<max;i++)
	*value++ = input->values[i];
  output->scale = input->scale;
  output->shift = input->shift;
  output->firstp = min;
  output->lastp = max-1;
  output->param = input->param;
}



/*--------------------------------------------------------------------------*/
/*
 * append a wrd into book
 */
/*--------------------------------------------------------------------------*/
void BookAppend(book, wrd)
     BOOK book;
     WRD wrd;
{
/* checking the inputs */

if ( book == NULL || wrd == NULL )
	perror("BookAppend: NULL input");	


  if(book->last == NULL) {
    book->first = wrd;
    book->last = book->first;
  }
  else {
    book->last->next = wrd;
    book->last = book->last->next;
  }
  book->size = book->size + 1;
  book->energy += wrd->coeff*wrd->coeff;
}


//----------------

void error(str)
char *str;
{
perror(str);
}

void error_option(str)
char *str;
{
perror("Option error\n");
}

void warning(str)
char *str;
{
printf(str);
}

