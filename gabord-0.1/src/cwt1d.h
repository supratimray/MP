/*..........................................................................*/
/*                                                                          */
/*       ------------------------Author Z. Zhang---------------------       */
/*         -------- (C) 1993 Copyright, All Right Reserved.--------         */
/*                                                                          */
/*..........................................................................*/
/****************************************************************************/
/*                                                                          */
/*  cwt1d.h    	 	Basic include file for the cwt files.    	    */
/*                                                                          */
/****************************************************************************/


/*
* transforms
*/
/* Changed from SIGNAL (*transform)[] to *transform[]  GD 7/11/93 */
extern SIGNAL *transform[MAX_NUM_SB];
/*
* cur and old_cur transform
*/
/*extern SIGNAL *cur_transform, *old_cur_transform; */

/* Replace cur_transform with a reference to the transform array */

#define cur_transform (transform[Current_Book])

/* temporary SIGNAL */
extern SIGNAL temporary;
/* 
 * number of transformations
 */
extern int TransAlloc[MAX_NUM_SB];
/* 
 * number of transformation cur and old_cur
 */
/*extern int *cur_TransAlloc;*/
#define cur_TransAlloc (TransAlloc[Current_Book])
extern int old_cur_TransAlloc;
/*
 * misc
 */
extern int ln2_dilatation;
extern int ln2_subsampling;
extern int octave_decomp; /* level of decomposition */
extern int SignalSize; /* real size of the initial signal */


extern int noct; /* noctave */
extern int nvoice; /* nvoice */
extern int O; /* number of band on the fine scales */
extern int F; /* number of fine scales */

extern double sigm1; /* filter parameter */
extern double sigm2; 
extern double sigm3;
extern int M; /* filter order */
extern int printing;
/*
 * end of cwt1d.h
 */

