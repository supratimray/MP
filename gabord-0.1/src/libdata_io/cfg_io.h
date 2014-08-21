/********************************************************************/
/*  Procedures for config files  input and output                   */
/* (C) 2001 Copyright Johns Hopkins University, All Right Reserved. */
/*  Piotr Franaszczuk JHU 2005                                      */

/* for ver. 1.3 5 mar  2002											*/
/* added intervals  10 Jan 2003										*/
/* added int_list  April 2004										*/
/********************************************************************/


#ifndef _CFG_IO_H_

#define _CFG_IO_H_
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <limits.h>

#define MAX_LINE 500
#define MAX_FILE_NAME 200
#define NAME_LEN 50

#ifdef WINDOWS
	#include <Windows.h>
#else
	#define BYTE unsigned char
#endif

// aliases for lists of double numbers
#define numbers_in  calib_in
#define numbers_out calib_out
#define double_list_in  calib_in
#define double_list_out calib_out


struct tag{
  char *name;
  void *addr;
  int (*fin)(char* str,void *buf);
  int (*fout)(void *buf, char* str);
};
struct section{
  char *name;
  struct tag* tags;
  int ntags;
  struct section *next;
};


#define set_lab_no(labels,n)  *(int*)labels=n
#define get_lab_no(labels)    *(int*)labels

#ifdef CPP
#define _DECL extern "C"
#else
#define _DECL extern
#endif

enum error_flag {INFO,WARNING,PERROR,ERROR,FERROR,DEBUG};

_DECL void print_error(enum error_flag f , char* pname, char *text,char*fname);
_DECL void mem_error(char *p,char* s);
_DECL int* iconv_numb(double x, double *buf);
_DECL double* conv_numb(double x, int *buf);
_DECL int* iconv_interv(double x, double *buf);
_DECL double* conv_interv(double x, int *buf);
_DECL int icomp_interv(int x, int *buf);
_DECL int comp_interv(double x, double *buf);
_DECL int intervals_in(char *s, void *buf);
_DECL int intervals_out(void* buf, char *s);
_DECL int float_in(char *s, void *buf);
_DECL int float_out(void* buf, char *s);
_DECL int int_in(char *s, void *buf);
_DECL int int_out(void* buf, char *s);
_DECL int double_in(char *s, void *buf);
_DECL int double_out(void* buf, char *s);
_DECL int short_in(char *s, void *buf);
_DECL int short_out(void* buf, char *s);
_DECL int str_in(char *s, void *buf);
_DECL int str_out(void* buf, char *s);
_DECL int enum_out(void* buf, char *s);
_DECL int enum_in(char *s, void* buf);
_DECL int labels_out(void* buf, char *s);
_DECL int labels_in(char *s, void* buf);
_DECL int calib_out(void* buf, char *s);
_DECL int calib_in(char *s, void* buf);
_DECL int int_list_out(void* buf, char *s);
_DECL int int_list_in(char *s, void* buf);
_DECL int function_out(void* buf, char *s);
_DECL int function_in(char *s, void* buf);
_DECL char** cross_labels(char **s);

_DECL struct section * create_section(char *name,struct tag *tags, int ntags);
_DECL int add_sections(struct section *list, struct section *add);
_DECL struct section* lookup_sec(struct section *head,char* name);
_DECL struct tag* lookup_tag(struct tag *tags, int n, char* name);
_DECL int write_tag(struct tag *t);
_DECL int write_sec(struct section *s);
_DECL int write_all(FILE* fil, struct section *head);
_DECL int skip_sec(char *name);   
_DECL int read_sec(struct  section *s); 
_DECL int read_all(FILE* fil, struct section *head,int cont);
_DECL int read_hdr(char*f, struct section *head);
_DECL int write_hdr(char*fi,char *fo, struct section *head);




#endif
