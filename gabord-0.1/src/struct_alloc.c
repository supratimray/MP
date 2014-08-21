/*--------------------------------------------------------------------------*/
/*  MPP signal processing program. 			                                */
/* (C) 1993 Copyright New York University, All Rights Reserved.			    */
/*							                                                */
/*  ------------------------------------------------------------            */
/*  Francois Bergeaud, Mike Orszag.                                         */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*                                                                          */
/*  struct_alloc.c   Functions which deal with the memory allocation        */
/*                   of BOOK structure                                      */
/*                                                                          */
/*--------------------------------------------------------------------------*/

#include "mpp.h"

/*--------------------------------------------------------------------------*/
/*
 *
 * Allocate a BOOK and retrun it
 *
 */
/*--------------------------------------------------------------------------*/

BOOK AllocBook()
{
  BOOK book;
  void init_book();

  if(!(book = (BOOK) (malloc(sizeof(struct book)))))
    error("Mem. alloc for STRUCTURE BOOK failed\n");

  init_book(book);

  return(book);
}

/*--------------------------------------------------------------------------*/
/*
 * initialize the book
 */
/*--------------------------------------------------------------------------*/
void init_book(book)
BOOK book;
{
   if (book == (BOOK)NULL)
	error("init_book(): null argument!");

   book->size = 0;
   book->energy =  0.0;
   book->sigen =  0.0;
   book->type = (int)NULL;
   book->smin = 0.;
   book->smax = 0.;
   book->sig_size= 0; 
   book->first = (WORD)NULL;
   book->last = (WORD)NULL;
}

/*--------------------------------------------------------------------------*/
/*
 * allocate WORD
 */
/*--------------------------------------------------------------------------*/
WORD AllocWord()
{
   WORD str;
   void init_word();

   if(!(str = (WORD) (malloc(sizeof(struct word)))))
	error("Mem. alloc for WORD failed\n");

   str->index = (INDEX)NULL;
   init_word(str);

   return(str);
}

/*--------------------------------------------------------------------------*/
/*
 * initialize the word
 */
/*--------------------------------------------------------------------------*/
void init_word(strt)
WORD strt;
{
    INDEX AllocIndex();
    void init_index();

    if (strt != (WORD)NULL)
	{
	strt->coeff = 0.0;
	strt->value = 0.0;
	strt->coeff1= 0.0;
	strt->value1= 0.0;
	strt->next = (WORD)NULL; /* wen */
	if (strt->index == (INDEX)NULL)
	   strt->index = AllocIndex();
	else
	   init_index(strt->index);
	strt->status = KEPT;
	}
    else
	warning("null arugment for init_word()!");
}

/*--------------------------------------------------------------------------*/
/*
 *  Allocation of index structure
 */
/*--------------------------------------------------------------------------*/
INDEX AllocIndex()
{
   INDEX indx;
   void init_index();

   if(!(indx = (INDEX) (malloc(sizeof(struct index)))))
	error("Mem. alloc for INDEX failed!");
   init_index(indx);
   return(indx);
}

/*--------------------------------------------------------------------------*/
/*
 *   Disallocation of index structure
 */
/*--------------------------------------------------------------------------*/
void init_index(index)
INDEX index;
{
   if (index != (INDEX)NULL)
	{
	index->id = 0.0;
	index->octave = 0.0;
	index->position = 0.0;
	index->phase = 0.0;
	}
   else
	warning("null argument for init_index()!");
}

/*--------------------------------------------------------------------------*/
/*
 *
 * deallocate the words field and initialize the other fields of the
 * structure book
 *
 */
/*--------------------------------------------------------------------------*/
clear_book(book)
BOOK book;
{
    WORD WordListFree();

    if (book == (BOOK)NULL)
	warning("null argument for clear_book()!");

	/* Free all the words in the book */
    WordListFree (book->first);

    init_book(book);
}

/*--------------------------------------------------------------------------*/
/*
 * deallocate the WORD structure
 */
/*--------------------------------------------------------------------------*/
WORD WordFree(strt)
WORD strt;
{
   INDEX IndexFree();

   if (strt == (WORD)NULL)
	warning("null argument for WordFree()!");

   strt->index = IndexFree(strt->index);
   free((char *)strt);
   return((WORD)NULL);
}

/*--------------------------------------------------------------------------*/
/*
 * deallocate the WORD list
 */
/*--------------------------------------------------------------------------*/
WORD WordListFree(word)
WORD word;
{
   WORD WordFree();

   if (word != (WORD)NULL)
        word->next = WordListFree(word->next);
   else
        return(word);
 
   word = WordFree(word);
   return((WORD)NULL);
}

/*--------------------------------------------------------------------------*/
/*
 * deallocate INDEX structure
 */
/*--------------------------------------------------------------------------*/
INDEX IndexFree(index)
INDEX index;
{
    if (index == (INDEX)NULL)
	warning("null argument for IndexFree()!");

    free((char *)index);
    return((INDEX)NULL);
}

/*--------------------------------------------------------------------------*/
/*
 * allocate index without initialization
 */
/*--------------------------------------------------------------------------*/
INDEX IndexAlloc()
{
    INDEX index;

    if ((index = (INDEX)malloc(sizeof(struct index)))
	== (INDEX)NULL)
	error("IndexAlloc(): mem. alloc. failed!");

    return(index);
}

/*--------------------------------------------------------------------------*/
/*
 * end of struct_alloc.c
 */
/*--------------------------------------------------------------------------*/

