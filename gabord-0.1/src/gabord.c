/****************************************************************************/
/*  GABOR signal processing program.                                        */
/* (C) 2001 Copyright Johns Hopkins University, All Right Reserved.         */
/*                                                                          */
/*  Christophe JOUNY: Based on Mallat/Zhang sources for Matching Pursuit    */
/*                    Adapted for continuous processing of EEGs             */
/*                    Added new energy threshold criteria                   */
/*  Piotr FRANASZCZUK: I/O procedures                                       */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*  gabord.c         Matching Pursuit decomposition adapted for processing  */
/*                   multiple signals and long datasets                     */
/*                                                                          */
/****************************************************************************/



#define structh

#include "data_io.h"

#include "mpp.h"
#include "structpath.h"
#include "cwt1d.h"
#ifndef WINDOWS
	#include <sys/time.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#ifdef WINDOWS
	#include <WinSock2.h>
	#include <string.h>
	#include <io.h>
#else
	#include <sys/times.h>
	#include <sys/param.h>
	#include <unistd.h>
	#include <strings.h>
#endif


extern struct inpdef *in;
extern struct outdef *ou, out_dflt;
extern int no_of_inp;
extern int cur_inp;

/*--------------------------------------------------------------------------*/

double *values=(double *)NULL;	//values of signal (temp var.)

SIGNAL signal=(SIGNAL)NULL;  	//signal
SIGNAL sigtmp=(SIGNAL)NULL;		//signal temp for decomp.
SIGNAL sigtmpi=(SIGNAL)NULL;	//    idem

double *cur_norm; // = (double *)NULL;

int cur_MaxOctave; /* Max octave used in the decomposition over gabor functions */
int cur_MinOctave; /* Min octave used in the decomposition over gabor functions */

double *pfG=(double *)NULL, *pfCE1=(double *)NULL;
double *pfCE2=(double *)NULL, *pfCE3=(double *)NULL, *pfC=(double *)NULL;
double *pfB = (double *)NULL;
int *pnAep=(int *)NULL, *pnIep=(int *)NULL;

int nBoundG=4; /* bound for the Gaussian */
int cur_l, cur_h; /* oversubsampling octaves for the fine grid */

/* input/output */
FILE *foutput;

/* library of books */
BOOK library[MAX_NUM_SB];
/* current book */
//BOOK cur_book;
/* previous book */
BOOK old_cur_book;
/* array of signals */
SIGNAL signals[MAX_NUM_SIGNAL];
/* current signal */
SIGNAL cur_signal;
/* current signal size */
int cur_sig_size=0;
/* filters */
SIGNAL * filter[MAX_NUM_SB];
/* current filter type */
int filter_type[MAX_NUM_SB];
/* shift octave and subsample octave */
int cur_shift_octave=0;
int cur_SOT=1;
int cur_SOF=1;
/* previous book */ 
SIGNAL *old_cur_filter;
/* number of filter */
int num_filter[MAX_NUM_SB];
/* previous number of filter */ 
int old_cur_num_filter = 0;
/*  global variable for transformation */
SIGNAL *(transform[MAX_NUM_SB]);   /* GD  7/11/93 */

int Current_Book = 0;
int Old_Book = 0;
/* number of transformation */
int TransAlloc[MAX_NUM_SB];
/* number of previous Transormation */
int old_cur_TransAlloc = 0;
/*  temporary Signal */
SIGNAL temporary;
/*  level of decomposition */
int octave_decomp = 0;
/*  size of the courant signal */
int SignalSize = 0;
/*  Constants defining the subsampling and dilatation rates */
int ln2_subsampling = 0;
int ln2_dilatation = 1;

double epslonG=1.0e-15; /* threshold for the Gaussian fuctions */

//Versioning number
const char version[]="1.4";

/**************************************************************************************/
/**************************************************************************************/


int main(argc, argv)
     int argc;
     char **argv;
	{
    int ShiftOctave=0;
    int SigSize;							//size of the signal analyzed
    int LnSigSize;							// power of 2 = to this size
    int SubsampleOctaveTime=2;
    int SubsampleOctaveFreq = 2;
    int SOT_build;
    int SOF_build;
    double sigma;
    double epslon;
    int old_cur_MinOctave = 0, old_cur_MaxOctave = 0; /* previous values of the levels of decomposition cur_MinOctave and cur_MaxOctave */

    double th_pct;						// criterion for terminate analysis based on prct of NRJ 
    int max_num_iter, stop_coh;						//  idem by number of atoms
    double th_atom_nrj, last_atom_nrj;                           // threshold of atoms energy and last atoms added nrj 
    double density=0.0;
    double energy=0.0;
    int cpct=0, citer=0, cth=0;
    ///////
    //extern FILE *foutput;
    char rbuffer[50];
    int sig_size = 32768;   			// max size of SIGNAL allowed
    int lndata, half_lndata, Lnlndata;	// length of data required by user, half of that and power of two
    //int blocksize;						// length of data block
    
    int sigN;							// signal size
	int nbinp;

    double *data;

    BOOK book;							//temp book for output
    WRD wrd;							//temp wrd from the book
    INDEX indx;							//temp index from the wrd
	
    int c, nb_chan, nbwind=0, windcount=0;
    int ichan;
    int n, nbvar;

	//DECLARATIONS FROM GBUILDBOOK

	/* the information for the fine grid */
	int curMaxLH; /* max of cur_l and cur_h */

	double SigEng, LamdaNoise; //factor=1.0, 
	long i, j, offset=0; //, ichan;
	int iter;
	static int nL=7, maxlh=0;
	double res_n, res_n1, orgN=-1.0;
	int l, h, nbc; 
	unsigned long flag=N_FLAG;
	//int FS = 200; //sampling frequency
  	int sb_index;
	short int *ds;
	float outval;
	double *df;


	union {double x; int p[2];} uniontemp;

	//Extern variables
	extern int no_of_out;
	

	// EXTERNAL FUNCTIONS

 	extern FILTER *GaborBuildFilter(), *GaborFreeFilter();
	extern double log2();
	extern SIGNAL *GaborDecomp();
	extern void change_signal();
	extern void sig_add_num();
	extern double sig_mean();

	/* building of the book */
	extern double *GaborGetGaussianArray();
	extern double *GaborGetCE1(), *GaborGetCE2(), *GaborGetCE3();
	extern double *GaborGetCArray(), *GaborGetB();
	extern int GaborGetBound();
	extern void GaborGetNAep(), GaborUpdateTrans();
	extern double farray_L2_sq_norm();
	extern BOOK AllocBook();
	extern BOOK init_book();
	extern WRD WordListFree();
	extern void GaborBuildBook();
	extern void GaborBuildBookOld();
	
	extern int find2power();
	extern void farray_copy();

//////// Section Definition

	struct tag xgbTag[]={
	  {"Percentage", &th_pct, double_in, double_out},
	  {"Energy", &th_atom_nrj, double_in, double_out},
	  {"Max_Iterations", &max_num_iter, int_in, int_out},
	  {"Coherence", &stop_coh, int_in, int_out}};

	struct section xgbSect={"GABOR_DECOMPOSITION", xgbTag, sizeof(xgbTag)/sizeof(struct tag), NULL};

////////////// HELP

	/*if (argc!=1)
	  {
	    printf("\n");
	    printf("Usage : gabord control_file\n");
	    printf("\n");
	    printf("\t-> control_file should be your control file with full path if necessary\n\n");
	    exit(0);
	  }*/

//////////// DEFAULTS SETTINGS

	th_pct =-100.0;
	max_num_iter = 10000;
	lndata = 32768;
	th_atom_nrj= -1;
	stop_coh=0;
	out_dflt.no_of_vars=-(3*max_num_iter+2);
	out_dflt.title=strdup("Book");

	foutput = stdout; //Redefined for compatibility with others files

/////////// DATA_IO INIT
//fprintf(stdout, "%i %i %i %i - %f %f %i\n", cpct, cth, stop_coh, citer, th_pct, th_atom_nrj,max_num_iter); 

	init_io(argc, argv, &xgbSect, 1);
	  
	// in->npoints=lndata;
	// in->win_len=lndata/in->freq;
	// in->win_size=lndata*nb_chan;

	//fprintf(stdout, "%i %i %i %i - %f %f %i\n", cpct, cth, stop_coh, citer, th_pct, th_atom_nrj,max_num_iter); 

	if (th_pct>-1.0) {cpct=1;}
	if (max_num_iter!=50000) {citer=1;}
	if (th_atom_nrj>-1.0) {cth=1;}

	if (th_pct<-1.0 & stop_coh==0 & th_atom_nrj<0 & citer==0) {citer=1; max_num_iter=1000;}

	ou[0].no_of_vars=max_num_iter*3+2;
	ou[0].var_len=1;
	
//fprintf(stdout, "%i %i %i %i - %f %f %i\n", cpct, cth, stop_coh, citer, th_pct, th_atom_nrj,max_num_iter); 

/* initializations     from int_loop.c  */

	for ( sb_index = 0; sb_index < MAX_NUM_SB; sb_index++) {
		filter[sb_index] = (SIGNAL *) NULL;
		num_filter[sb_index] = 0;
		TransAlloc[sb_index] = 0;
		filter_type[sb_index] = -1;
		}
/* allocation of book */
    for (i=0;i<MAX_NUM_SB;i++)
		{
        library[i] = AllocBook();
        library[i]->id=i;
		}
/* allocation of signals  */
	for (i=0;i<MAX_NUM_SIGNAL;i++)
        signals[i] = new_struct_signal();
	
/*cur_book = library[0];*/
	old_cur_book = cur_book;

/* assign the current signal */
	cur_signal = signals[0];


	temporary = new_struct_signal();
	old_cur_filter = filter[0];
	Current_Book = 0;

// Others init

	ds = (short int *)malloc(sizeof(short int));
	df = (double *)malloc(sizeof(double));

////

	signal = new_signal(sig_size);

	Lnlndata = find2power(lndata);
	lndata = 1<<Lnlndata;

/// START LOOPING

	for (nbinp=0; nbinp<no_of_inp; nbinp++)
		{
		cur_inp=nbinp;
		
		nb_chan=in[nbinp].no_of_chans;
		lndata=in[nbinp].npoints;
	
		if (lndata>32768)
	  	{
	//	char rbuffer[50];
	  	lndata=32768;
		sprintf(rbuffer, "Window Length Changed to %i points\n", lndata);
	  	print_error(WARNING, "main", rbuffer, NULL);
	  	
	
		new_inp_size(in+nbinp, lndata);
		}


		while (n=read_next_win(in+nbinp, -1, 0))
			{

			windcount++;
		//fprintf(stderr, "%i - %lf\n", windcount-1, density);


		for (ichan=0; ichan<nb_chan; ichan++)
			{

			signal->size_alloca = sig_size;
			signal->size = lndata;
			signal->shift = 0;
			signal->scale = 1;
			signal->firstp = 0;
			signal->lastp= signal->size - 1;
			signal->param = 1;

			//Read Channel Data
			nbc=get_double_chan(in+nbinp, signal->values, n, ichan);

  			sig_add_num(signal, -1*sig_mean(signal));		// substract mean to center around zero

			orgN = (double) signal->size;
			sigN = signal->size;
			LnSigSize = find2power(signal->size);
			cur_MaxOctave = LnSigSize - 1;
			cur_SOT=1; cur_SOF=1;
			ShiftOctave = 0;

			/* previous values for the decomposition */

	    		old_cur_MinOctave = cur_MinOctave; 
    			old_cur_MaxOctave = cur_MaxOctave;
		    	cur_MinOctave = 1;
    			SOT_build = cur_SOT + 2;
    			SOF_build = cur_SOF + 2;

		    	SigSize = 1<<LnSigSize;
    			if (SigSize != sigN)         				       // signal resize to a power of 2
				{
   				values=(double *)malloc(sizeof(double)*sigN);
		    		if (values  ==(double *)NULL)
					perror("GDecomp(): mem. alloc. failed!");
				farray_copy(signal->values,sigN,values);
				change_signal(signal,SigSize);
				for (i=0;i<sigN;i++)
		    			signal->values[i] = values[i];
				for (i=SigSize;i<sigN;i++)
	    				signal->values[i] = 0.0;
				free((char *)values);
				sigN = SigSize;
				}

			
		    	if (SubsampleOctaveTime > LnSigSize)
				SubsampleOctaveTime = LnSigSize;
		    	if (SubsampleOctaveFreq> LnSigSize)
				SubsampleOctaveFreq = LnSigSize;
		    	if (ShiftOctave == 0)
				ShiftOctave = (LnSigSize+SubsampleOctaveTime-SubsampleOctaveFreq+1)/2;
	
			 sigma = (double)(1.0/sqrt((double)(2.0*M_PI))); 	/* 1/sqrt(2*pi) */

    			if (cur_num_filter!=LnSigSize+1 || cur_filter_type != NEWGABOR)
				{										/* free the gabor filter */
				if (cur_filter != (FILTER *)NULL)
					cur_filter = GaborFreeFilter(cur_filter,cur_num_filter);
				if (cur_transform != (SIGNAL *)NULL)
					cur_transform = GaborFreeFilter(cur_transform, cur_TransAlloc);
				}

		    	if ( ( old_cur_MinOctave != cur_MinOctave ) || ( old_cur_MaxOctave != cur_MaxOctave ) ) /* octave have changed */
				{
				if (cur_transform != (SIGNAL *)NULL)
					cur_transform = GaborFreeFilter(cur_transform, cur_TransAlloc);
				}
	
    			if (cur_filter == (FILTER *)NULL)
				{
				cur_filter = GaborBuildFilter( LnSigSize, SigSize, sigma,&cur_norm);
				cur_num_filter = LnSigSize+1;
				cur_filter_type = NEWGABOR;
				}

													/* Do the decomposition */
		    	if (sigtmp == (SIGNAL)NULL)
				sigtmp = new_signal(SigSize<<1);
		    	else if (sigtmp->size != SigSize)
				change_signal(sigtmp,SigSize<<1);

		    	for (i=0;i<SigSize;i++)
				sigtmp->values[i] = signal->values[i];
				
		    	cur_transform = GaborDecomp(cur_transform,sigtmp,cur_filter,
					SubsampleOctaveTime,SubsampleOctaveFreq,
					cur_MinOctave, cur_MaxOctave,ShiftOctave);
		/*
 		* set the last ShiftOctave and SubsampleOctaveTime and SigSize
 		* to avoid recomputation of the filters
 		*/
			    cur_TransAlloc = cur_MaxOctave - cur_MinOctave +3;
			    cur_shift_octave = ShiftOctave;
			    cur_SOT = SubsampleOctaveTime;
			    cur_SOF = SubsampleOctaveFreq;
			    cur_sig_size = SigSize;
			    cur_signal = signal;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

		    epslon = epslonG;
		    l=cur_SOT+2;
		    h=cur_SOF+2;
		
		    cur_l = l-1;
		    cur_h = h-1;
		
		    if (cur_transform == (SIGNAL *)NULL || cur_filter_type != NEWGABOR)
				perror("Run gdecomp first!");
		 
		    if (cur_book==(BOOK)NULL)
		        cur_book = AllocBook();
		    if (cur_book->first != (WRD)NULL)
		        cur_book->first = WordListFree(cur_book->first);
		    init_book(cur_book);
 
		    if (cur_transform[0] == (SIGNAL)NULL)
				perror("Internal Error!");

    		SigEng = farray_L2_sq_norm(cur_transform[0]->values, cur_transform[0]->size>>1);

		    cur_book->sigen = res_n = SigEng;
		    cur_book->sig_size = sigN;
		    cur_book->type = NEWGABOR;

		    if (find2power(cur_sig_size)!=nL)
				{
				nL = find2power(cur_sig_size);
				if (pnAep!=(int *)NULL)
				    {
				    free((char *)pnAep);
				    pnAep = (int *)NULL;
				    }
				if (pnIep!=(int *)NULL)
				    {
				    free((char *)pnIep);
				    pnIep = (int *)NULL;
				    }
				if (pfG!=(double *)NULL)
				    {
				    free((char *)pfG);
				    pfG=(double *)NULL;
				    }
				if (pfCE1!=(double *)NULL)
				    {
				    free((char *)pfCE1);
				    pfCE1=(double *)NULL;
				    }
				if (pfC!=(double *)NULL)
				    {
				    free((char *)pfC);
				    pfC=(double *)NULL;
				    }
				if (pfB!=(double *)NULL)
				    {
				    free((char *)pfB);
				    pfB = (double *)NULL;
				    }
				}
		
			if ((flag&N_FLAG)==N_FLAG)
		    	{
		    	if (maxlh!=MAX(MAX(cur_SOT-1,cur_SOF-1),MAX(l-1,h-1)))
				{
				curMaxLH = maxlh = MAX(MAX(cur_SOT-1,cur_SOF-1),MAX(l-1,h-1));
				if (pfG != (double *)NULL)
			    	{
			    	free((char *)pfG);
			    	pfG = (double *)NULL;
			    	}
				}
    
		    if (pfB == (double *)NULL)
				pfB = GaborGetB(nL,epslon);
		
		    if (pfC==(double *)NULL)
				pfC = GaborGetCArray(nL);
		
		    if (pfG==(double *)NULL)
				{
				if (pnIep==(int *)NULL)
				    pnIep = (int *)malloc(sizeof(int)*nL);
				if (pnAep==(int *)NULL)
				    pnAep = (int *)malloc(sizeof(int)*nL);
				if (pnAep==(int *)NULL || pnIep==(int *)NULL)
				    perror("mem. alloc. failed!");
				GaborGetNAep(nL,maxlh,epslon,pnAep,pnIep);
				pfG = GaborGetGaussianArray(nL,pnAep,pnIep[nL-1],maxlh);
				}
		
		    if (pfCE1==(double *)NULL)
				pfCE1 = GaborGetCE1(nL);
		    if (pfCE2==(double *)NULL)
				pfCE2 = GaborGetCE2(nL);
		    if (pfCE3==(double *)NULL)
				pfCE3 = GaborGetCE3(nL);
		    }

///////////////////////////////////////////////////////////    STOP CRITERION

/* compute the lamda square for white noise */

		    if (SigEng < 0.0)
				perror("GaborBuildBook(): empty signal!");
		    if (cur_l != 3 || cur_h != 3)
				LamdaNoise = 0.0; /* this lambda is not available yet */
		    else
				LamdaNoise = (2.076747-0.089091*log(orgN))*sqrt(log(orgN)/orgN);


		    res_n1 = 0.0;
		    last_atom_nrj = SigEng;             // to pass the first test : add xtof

		    for (iter=0;iter<max_num_iter;iter++)			// Iteration criteria
		      {
			
			if (cpct)
			  if (cur_book->energy>=SigEng*th_pct/100.0)	                // Threshold criteria
			      break;
				
			if (stop_coh)
			  if (sqrt(1.0-res_n1/res_n)<LamdaNoise)		// Coherence criteria
			      break;
			
			if (cth)
			  {
//			    fprintf(stdout, "NRJ = %lf\n", last_atom_nrj);
			  if (last_atom_nrj < th_atom_nrj)			// Atoms NRJ criteria
			      break;
			  }
			GaborBuildBook(cur_transform,cur_filter,cur_book, 
							cur_MinOctave,
							cur_MaxOctave, 
							nL,
							cur_SOT,cur_SOF,l,h,
		                				pnIep,pfC,pfB,pfG,pfCE1,
		                				cur_norm);

			last_atom_nrj = cur_book->last->coeff*cur_book->last->coeff;

				
	//		if ((flag&D_FLAG)==D_FLAG)
		//	  {
			    if (iter>0)
			      res_n = res_n1;
			    res_n1 = res_n - last_atom_nrj;
			//  }

		      }

		book = cur_book;

//////////////WRITING THE BOOK
		
		nbvar=3*book->size+2;
			  
		uniontemp.p[0]=windcount;
		uniontemp.p[1]=book->size;
		ou[0].write_chan_data(&ou[0],  &uniontemp.x, 1, ichan, nbvar);

	  	energy  = sqrt(book->sigen); //should not be sqrt... stay for now for backward coherence in book files
	  	ou[0].write_chan_data(&ou[0], &energy, 1, ichan, 0);  // Signal Energy

		wrd = book->first;
		i = 0;

		while (wrd != NULL)
		  {
		    union {double x; short int ps[4];} uniontemp2;

		    indx = wrd->index;
		  
		    uniontemp2.ps[0]=i++;
		    uniontemp2.ps[1]=indx->octave;
		    uniontemp2.ps[2]=indx->id;
		    uniontemp2.ps[3]=indx->position;	 	    

		    ou[0].write_chan_data(&ou[0], &uniontemp2.x, 1, ichan, 0);
		    ou[0].write_chan_data(&ou[0], &(wrd->coeff), 1, ichan, 0);
		    ou[0].write_chan_data(&ou[0], &(indx->phase), 1, ichan, 0);
	      			
		    wrd = wrd->next;
		  }

			  
			  
////////////WRITING GAD FILE
		  
		if (no_of_out==2)
			{
		
			density = 2.0 * ((double)book->size) / book->sig_size;
	
			ou[1].write_chan_data(&ou[1], &density, 1, ichan, 4); //Density
	
			energy  = sqrt(book->energy);
			ou[1].write_chan_data(&ou[1], &energy, 1, ichan, 0);  // Book Energy
			energy  = sqrt(book->sigen);
			ou[1].write_chan_data(&ou[1], &energy, 1, ichan, 0);  // Signal Energy
	
			uniontemp.p[0]=windcount;
			uniontemp.p[1]=book->size;
			ou[1].write_chan_data(&ou[1],  &uniontemp.x, 1, ichan, 0);  // Counter & book size

			}


			}

		}
	}
	if(windcount<1)print_error(ERROR,"main","Nothing computed",NULL);
	close_io();
	//	print_error(INFO, "main", "Done", NULL);
	exit(0);
	}
	

	

/*--------------------------------------------------------------------------*/
 /*
 * end of gabord.c
 */
/*--------------------------------------------------------------------------*/
