
/********************************************************************/
/*  Procedures for config files  input and output                   */
/* (C) 2001 Copyright Johns Hopkins University, All Right Reserved. */
/*  Piotr Franaszczuk JHU 2001,2002                                 */ 
/*  ver. 1.5              17 jul  2002                              */
/*  Jan 2003  added interv functions and alias for calib            */
/* 	$Id: cfg_io.c,v 1.5 2007/05/03 20:58:32 pfranasz Exp $	        */
/********************************************************************/
#ifdef WINDOWS
	#include <stddef.h>
	#include <string.h>
	#include <malloc.h>
char * strndup (const char *s, size_t n)
{
  char *result;
  size_t len = strlen (s);

  if (n < len)
    len = n;

  result = (char *) malloc (len + 1);
  if (!result)
    return 0;

  result[len] = '\0';
  return (char *) memcpy (result, s, len);
}
#endif

#ifndef lint
static char vcid[] = "$Id: cfg_io.c,v 1.5 2007/05/03 20:58:32 pfranasz Exp $";
#endif /* lint */
#define _GNU_SOURCE
#include "cfg_io.h"

const char cfg_ver[]="1.5";
extern char *progname,*host;
#ifndef WINDOWS
	extern int errno;
#endif
char *del_list[100];

int del_no=0;

char* s_fmt[]= {"Little_Endian","Big_Endian",
	        "text",
		"char","byte",
		"int12",
		"short","ushort","edf",
		"int","uint","float",
		"double",
		"serial","parallel",
		"selections","files"
};
#define S_FMT_SIZE sizeof(s_fmt)/sizeof(char*)

#ifdef _WIN32
#define rint(x) floor(x+0.5)
#endif



char section_char='%';
char *delimiters="=,:";    /* has to be global */
FILE *inp_ops_fil, *out_ops_fil;
char buf_line[MAX_LINE];
int copy_ops=0;            /* default no copy in read_all */
int print_warning=1;

void printerror_(char * f,char * p, char * t, char *fname,
			int lf, int lp, int lt, int lfn)
{
	// F77 wrapper it has to be consistent with enum error_flag in cfg_io.h
	static char flags[6]={'I','W','P','E','F','D'};
	int i; 
	char* pp,*tt,*ff;
	ff=strndup(fname,lfn);
	pp=strndup(p,lp);
	tt=strndup(t,lt);
	for(i=0;i<6;i++)
	    if(*f==flags[i]) break;
	print_error(i,pp,tt,ff);
	free(tt);
	free(pp);
	free(ff);
	return;
}

void print_error(enum error_flag f, char* procnam, char * text, char* fname)
{
  char buf[1000];
  switch(f)
    {
    case INFO: 
      fprintf(stderr,"%s:%s\t%s: %s\n",host,progname,procnam,text);
      return;
    case DEBUG:
#ifdef DEB
      fprintf(stderr,"%s:%s\tDEBUG:%s: %s\n",host,progname,procnam,text);
      ffflush(stderr);
#endif
      return;
    case WARNING:
      if(!print_warning)return;
      fprintf(stderr,"%s:%s\tWARNING:%s: %s\n",host,progname,procnam,text);
      fflush(stderr);
      return;
    case ERROR:
       fprintf(stderr,"%s:%s\tERROR:%s: %s\n",host,progname,procnam,text);
       break;
    case PERROR:
       sprintf(buf,"%s:%s\tERROR:%s: %s",host,progname,procnam,text);
       perror(buf);
       break;
    case FERROR:
	fprintf(stderr,"%s:%s\tERROR:%s: %s file: %s\n",host,progname,procnam,text,fname);
	if(errno)
	  perror(NULL);
      break;
    default:
	 fprintf(stderr,"%s:%s\tERROR:UNKNOWN:%s: %s\n",host,progname,procnam,text);
    } 
  exit(3);
}

void mem_error(char *p,char* s)
{ char erbuf[200];
 
  sprintf(erbuf,"memory allocation of %s",s);
  print_error(ERROR,p,erbuf,NULL);
 
}
int * iconv_interv(double x,double*interv)
{
  // converts interval list from double to int
  // i=interv*x 
  int i,n,*ret;
  if(!interv)return NULL;
  n=(int)rint(interv[0]);
  if(n<1)return NULL;
  if(!(ret=calloc(2*n+1,sizeof(int))))mem_error("iconv_interv","ret");
  ret[0]=n;
  for(i=1;i<=2*n;i++)
    ret[i]=(int)rint(interv[i]*x);
  return ret;
}
 
double * conv_numb(double x, int *numb)
{
  // converts number list from int to double 
  // i=numb*x 
  int i,n;
  double *ret;
  if(!numb)return NULL;
  n=numb[0];
  if(n<1)return NULL;
  if(!(ret=calloc(n+1,sizeof(double))))mem_error("conv_numb"," ");
  ret[0]=n;
  for(i=1;i<=n;i++)
    ret[i]=numb[i]*x;
  return ret;
}

int * iconv_numb(double x,double*numb)
{
  // converts number list from double to int
  // i=numb*x 
  int i,n,*ret;
  if(!numb)return NULL;
  n=(int)rint(numb[0]);
  if(n<1)return NULL;
  if(!(ret=calloc(n+1,sizeof(int))))mem_error("iconv_numb"," ");
  ret[0]=n;
  for(i=1;i<=n;i++)
    ret[i]=(int)rint(numb[i]*x);
  return ret;
}
 
double * conv_interv(double x, int *interv)
{
  // converts interval list from int to double 
  // i=interv*x 
  int i,n;
  double *ret;
  if(!interv)return NULL;
  n=interv[0];
  if(n<1)return NULL;
  if(!(ret=calloc(2*n+1,sizeof(double))))mem_error("conv_interv"," ");
  ret[0]=n;
  for(i=1;i<=2*n;i++)
    ret[i]=interv[i]*x;
  return ret;
}
int icomp_interv(int x, int*interv)
{
  // compares if x is in any of intervals in interv (integer version)
  // returns i if it is in i-th interval; 0 if it is not in any
  // interv is stored as read by intervals_in
  // if interv==NULL returns -1

  int i,n;
  if(!interv)return -1;
  n=interv[0]*2;
  for(i=1;i<n;i+=2)
    {
      if(x<interv[i] || x> interv[i+1])continue;
      return (i+1)/2;
    }
  return 0;
}

int comp_interv(double x, double*interv)
{
  // compares if x is in any of intervals in interv
  // returns i if it is in i-th interval; 0 if it is not in any
  // interv is stored as read by intervals_in
  // if interv==NULL returns -1

  int i,n;
  if(!interv)return -1;
  n=(int)rint(interv[0])*2;
  for(i=1;i<n;i+=2)
    {
      if(x<interv[i] || x> interv[i+1])continue;
      return (i+1)/2;
    }
  return 0;
}
  
int intervals_in(char* s, void *buf)
{
  /* input for list of pairs of doubles defining intervals */
  /* it should NEVER be called directly, only from read_sec !!!  */
  /* buf is created here the size is taken from input         */
  /* buf[0] contains number of non-overlapping intervals (as int)   */

  int i,n,ret,k;
  double *bb;
  double a,b;
  
  if(!s)goto err;
   
  n=atoi(s);
  if(n<1)return 0;
  if((bb=(double *)malloc((2*n+1)*sizeof(double)))==NULL)mem_error("intervals_in","bb");
  bb++;
  s=NULL;
  k=0;
  for(i=0;i<n;i++)
    {   
      while(s==NULL)
	{
	  if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)goto err;
	  s=strtok(buf_line,")");
	}
     
      if(sscanf(s,"(%lf,%lf",&a,&b)!=2)goto err;

      if(a<b) // use only non-empty intervals
	{
	  bb[k++]=a;
	  bb[k++]=b;
	}
      s=strtok(NULL,")");	 	  
    }
  bb--;
  bb[0]=k/2;
  *(double**)buf=bb;
  return 0;
 err: 
  print_error(ERROR,"intervals_in","parsing input string",NULL);
  
    
}   
      

int intervals_out(void *buf, char*s)
{
  int i,n;
  double *bb=*(double**)buf;
  if(bb==NULL || (n=(int)rint(bb[0]))<1){
    strcpy(s,"0");
    return 1;}
  fprintf(out_ops_fil,"%d\n",n);
  n*=2;
  for(i=1;i<n-2;i+=2)
    fprintf(out_ops_fil,"(%.16g,%.16g)\n",bb[i],bb[i+1]);
  fprintf(out_ops_fil,"(%.16g,",bb[n-1]);
  sprintf(s,"%.16g)",bb[n]); /* last is written outside this proc */
  if(copy_ops)
    { char *b;
    // skip all 
    b=strtok(NULL,delimiters);
    n=atoi(b);
    for(i=0;i<n;i++)
      fscanf(inp_ops_fil,"%*[^\n]\n");
    }
  return strlen(s);
}  


int float_in(char *s, void *buf)
{
  *(float*)buf =(float)atof(s);
  return 0;
}
int float_out(void* buf, char *s)
{ return sprintf(s,"%.8g",*(float*)buf);}


int int_in(char *s, void *buf)
{
  *(int*)buf=atoi(s);
  return 0;
}

int int_out(void* buf, char *s)
{ return sprintf(s,"%i",*(int*)buf);}

int double_in(char *s, void *buf)
{
  *(double*)buf=atof(s);
  return 0;
}
int double_out(void* buf, char *s)
{ return sprintf(s,"%.16g",*(double*)buf);}


int short_in(char *s, void *buf)
{
  *(short*)buf=(short)atoi(s);
  return 0;
}

int short_out(void* buf, char *s)
{ return sprintf(s,"%hi",*(short*)buf);}


int str_in(char *s, void *buf)
{
 
  if(*(char**)buf)
      free(*(char**)buf);
  if(!(*(char**)buf=strdup(s)))mem_error("str_in","buf");
  
  return strlen(s);
}

int str_out(void* buf, char *s)
{ 
  if(*(char**)buf)
      strncpy(s,*(char**)buf,MAX_LINE);
  else
       strcpy(s,".");
     
  return strlen(s);
}

int enum_out(void* buf, char *s)
{
  int i=*(int*)buf;
  if(i<0 || i >=S_FMT_SIZE)
    strcpy(s,"UNKNOWN");
  else
    strcpy(s,s_fmt[i]);
  return strlen(s);
}

int enum_in(char *s, void* buf)
{
  int i;
  for(i=0;i<S_FMT_SIZE;i++)
    if(!strcmp(s,s_fmt[i]))break;
  if(i==S_FMT_SIZE)i=-1;
  *(int*)buf=i;
  return i;
}

int labels_in(char* s, void *buf)
{
  /* input for list of strings delimited by characters in global var delimiters or new lines */
  /* it should NEVER be called directly, only from read_sec !!!                           */
  /* buf is created here the size is taken from first                                     */

  int i,n,ret;
  char **bb;
  if(!s)
    {
      print_error(ERROR,"labels_in","Missing no of labels",NULL);
     
    }
  n=atoi(s);
  if(n<1)return 0;
  s=NULL;
  if((bb=(char**)calloc(n+1,sizeof(char*)))==NULL)mem_error("labels_in","bb");
  set_lab_no(bb,n);
  bb++;
  for(i=0;i<n;i++)
    {    
      while(s==NULL)
	{
	  if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)
	    { char erbuf[100];
	    sprintf(erbuf,"Reading label %i",i+1);
	     print_error(ERROR,"labels_in",erbuf,NULL);
	    }
	  s=strtok(buf_line,delimiters);
	}
      str_in(s,bb+i);
      s=strtok(NULL,delimiters);
    }
  *(char***)buf=bb-1;	
  return 0;
}   
      
int labels_out(void *buf, char*s)
{
  int i,n,k;
  char **bb=*(char***)buf;
  char *b;
  if(copy_ops && (b=strtok(NULL,delimiters)))
    { 
    // skip all labels
	k=atoi(b);
	for(i=0;i<k;i++)
	  fscanf(inp_ops_fil,"%*[^\n]\n");
    }
  if(bb==NULL ||(n=get_lab_no(bb))<1 ){
    strcpy(s,"0");
    return 1;}

  fprintf(out_ops_fil,"%d\n",n);
  for(i=1;i<n;i++)
    fprintf(out_ops_fil,"%s\n",bb[i]);
      
  if(!bb[n])mem_error("labels_out","bb[n]");

  
  strncpy(s,bb[n],MAX_LINE); /* last is written outside this proc */
  return strlen(s);
}
   
int calib_in(char* s, void *buf)
{
  /* input for list of doubles delimited by characters in global var delimiters or new lines */
  /* it should NEVER be called directly, only from read_sec !!!                           */
  /* buf is created here the size is taken from input                                     */

  int i,n,ret;
  double *bb;
  if(!s)
    {
       print_error(ERROR,"calib_in","Missing no of calib factors",NULL);
      
    }
  n=atoi(s);
  if(n<1)return 0;
  s=NULL;
  if((bb=(double *)malloc((n+1)*sizeof(double)))==NULL)mem_error("calib_in","bb");
  bb[0]=n;
  bb++;
  for(i=0;i<n;i++)
    {   
      while(s==NULL)
	{
	  if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)
	    { char erbuf[100];
	    sprintf(erbuf,"Reading item %i",i+1);
	     print_error(ERROR,"calib_in",erbuf,NULL);
	    };
	  s=strtok(buf_line,delimiters);
	}
      double_in(s,bb+i);
      s=strtok(NULL,delimiters);
    }
  *(double**)buf=bb-1;
  return 0;
}   
      
int calib_out(void *buf, char*s)
{
  int i,n;
  double *bb=*(double**)buf;
 if(copy_ops)
    { char *b;
    // skip all 
    b=strtok(NULL,delimiters);
    n=atoi(b);
    for(i=0;i<n;i++)
      fscanf(inp_ops_fil,"%*[^\n]\n");
    }
  if(bb==NULL || (n=(int)rint(bb[0]))<1){
    strcpy(s,"0");
    return 1;}
  fprintf(out_ops_fil,"%d\n",n);
  for(i=1;i<n;i++)
    fprintf(out_ops_fil,"%.16g\n",bb[i]);
  double_out(bb+n,s); /* last is written outside this proc */
  
  return strlen(s);
}  
int int_list_in(char* s, void *buf)
{
  /* input for list of integers delimited by characters in global var delimiters or new lines*/
  /* it should NEVER be called directly, only from read_sec !!!                           */
  /* buf is created here the size is taken from input                                     */

  int i,n,ret;
  int *bb;
  if(!s)
    {
       print_error(ERROR,"int_list_in","Missing no of integers in list",NULL);
      
    }
  n=atoi(s);
  if(n<1)return 0;
  s=NULL;
  if((bb=(int *)malloc((n+1)*sizeof(int)))==NULL)mem_error("int_list__in","bb");
  bb[0]=n;
  bb++;
  for(i=0;i<n;i++)
    {   
      while(s==NULL)
	{
	  if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)
	     { char erbuf[100];
	     sprintf(erbuf,"Reading item %i",i+1);
	     print_error(ERROR,"int_list_in",erbuf,NULL);
	    }
	  s=strtok(buf_line,delimiters);
	}
      int_in(s,bb+i);
      s=strtok(NULL,delimiters);
    }
  *(int**)buf=bb-1;
  return 0;
}   
      
int int_list_out(void *buf, char*s)
{
  int i,n;
  int *bb=*(int**)buf;
if(copy_ops)
    { char *b;
    // skip all 
    b=strtok(NULL,delimiters);
    n=atoi(b);
    for(i=0;i<n;i++)
      fscanf(inp_ops_fil,"%*[^\n]\n");
    }
  if(bb==NULL || (n=bb[0])<1){
    strcpy(s,"0");
    return 1;}
  fprintf(out_ops_fil,"%d\n",n);
  for(i=1;i<n;i++)
    fprintf(out_ops_fil,"%d\n",bb[i]);
  int_out(bb+n,s); /* last is written outside this proc */
  
  return strlen(s);
}  


struct section * create_section(char *name,struct tag *tags, int ntags)
{
  struct section *this;
  this=(struct section*)malloc(sizeof(struct section));
  if(!this)mem_error("create_section","this");
  this->name=name;
  this->tags=tags;
  this->ntags=ntags;
  this->next=NULL;
  return this;
}
int add_sections(struct section *list, struct section *add)
{
  // at end
  struct section *this=list;
  while(this->next) this=this->next;
  this->next=add;
  return 1;

}

struct section* lookup_sec(struct section *head,char* name)
{
  struct section *this=head;
  if(this->name==NULL)return this; // to read parameters without section name;
  while(this && strcmp(name,this->name)) this=this->next;
  return this;
}
int tag_count;

struct tag* lookup_tag(struct tag *tags, int n, char* name)
{
  struct tag *this=tags; 
  int i;
  for(i=0;i<n;i++,this++)
  {
    if(!strcmp(name,this->name))break;
    if(this->name[0]=='*' && i==tag_count)
	return this;
	
  }
  if(i<n)return this;
  else return NULL;
}

    


int write_tag(struct tag *t)
{ int ret;
 if(t->fout==NULL || t->addr==NULL )return 1;  // do not output this tag
 fprintf(out_ops_fil,"%s=",t->name); 
 if((ret=(*t->fout)(t->addr,buf_line))<1)return ret;
 return fprintf(out_ops_fil,"%s\n",buf_line);

}

int write_sec(struct section *s)
{  int ret,i;
 struct tag *t=s->tags;
 fprintf(out_ops_fil,"%%%s\n",s->name);
 for(i=0;i<s->ntags;i++,t++)
   if((ret=write_tag(t))<1)return ret;
 return i;
}

int write_all(FILE* fil, struct section *head)
{
  int ret,n=0;
  struct section *s=head;
  out_ops_fil=fil;
  while(s)
    {
      if((ret=write_sec(s))<1)return ret;
      n++;
      s=s->next;
    }
  
  return n;
}

int skip_sec(char *name)
{
  int ret,i,notdel=1;
  do{
    while(buf_line[0]!=section_char)
      {
	if(copy_ops && notdel )fprintf(out_ops_fil,"%s\n",buf_line);
	if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)return ret;
      }
    strcpy(name,buf_line+1);
    if(copy_ops && del_no>0)
      {
	notdel=1;
	for(i=0;i<del_no;i++)
	  {
	    if(!strcmp(name,del_list[i]))
	      {notdel=0;break;}
	  }
        if(notdel)break;
        if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)return ret;
      }
   
  }while(copy_ops && !notdel);
  if(copy_ops)fprintf(out_ops_fil,"%s\n",buf_line);
  return strlen(name);
}
  
int read_sec(struct  section *s)
{
  int ret,i,n=0;
  struct tag *t;
  char  *b;
  int(**fout)(void*,char*);

  if(s->name && (ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)return ret;
  
  
  fout=malloc(s->ntags*sizeof(int(*)(void*,char*)));
  if(!fout)mem_error("read_sec","fout");
  tag_count=0;
  for(i=0;i<s->ntags;i++)
	  fout[i]=s->tags[i].fout;
  
  while(buf_line[0]!=section_char)
    {
      if(buf_line[0]!='#')
	{ char * buf=strdup(buf_line);
	  b=strtok(buf_line,delimiters);
	  if(b)
	    {
	      if((t=lookup_tag(s->tags,s->ntags,b)))
		{/* tag in section */
		 	
		  if(copy_ops)
		    {
			write_tag(t);
			fflush(out_ops_fil);
			t->fout=NULL;  // to avoid writing it twice
		    }
		  else
		    {
		       if(t->name[0]!='*')// not wild card tag
			 b=strtok(NULL,delimiters);
		       else
			 tag_count++;	
		      if(!b)
			{char erbuf[100];
			 sprintf(erbuf,"Empty tag %s in section %s",buf_line,s->name);
			 print_error(ERROR,"read_sec",erbuf,NULL);
		       }
		      if(t->fin && t->addr)(*t->fin)(b,t->addr);
			 
		    }
		  n++;	  
		}
	      else
		  {
		      // not recognized tag
		     if(!copy_ops)
		       {char erbuf[100];
			 sprintf(erbuf,"tag %s in section %s not processed",b,s->name);
			 print_error(WARNING,"read_sec",erbuf,NULL);
		       }
		     else 
		        {fprintf(out_ops_fil,"%s\n",buf);fflush(out_ops_fil);}

		  }
	    }
	  if(buf)free(buf);
	}
      else
	if(copy_ops)
	    {fprintf(out_ops_fil,"%s\n",buf_line);fflush(out_ops_fil);}

      if((ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line))!=1)break;
    }
  if(n!=s->ntags)
      {
      if(copy_ops)
	  {
	     char *buf =strdup(buf_line);
	     for(i=0;i<s->ntags;i++)
		 {  
		    write_tag(s->tags+i);
		    fflush(out_ops_fil);
		 }
	     strcpy(buf_line,buf);
	     free(buf);
	  }
      else
	{char erbuf[100];
	      sprintf(erbuf,"not all tags defined in section %s",s->name);
	      print_error(WARNING,"read_sec",erbuf,NULL);
		       }
	   
     }
 // copy fout back
  for(i=0;i<s->ntags;i++)
	  s->tags[i].fout=fout[i];  
  free(fout);
  return n+1;
}
   
int read_all(FILE* fil, struct section *head, int cont)
{
  struct section *s;
  int ret=1,n=0;
  char name[MAX_LINE];
  inp_ops_fil=fil;
  if(!cont)ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line);// for cont=1 buf_line already contains sec name
  while(!feof(fil))
    {
      if(ret<1)return ret;
      if(head->name!=NULL && (ret=skip_sec(name))<1)return n;
      if((s=lookup_sec(head,name)))
	{
	  ret=read_sec(s);
	  n++;
	}
      else
	  {
	      if(!copy_ops)
		{ char erbuf[100];
		  sprintf(erbuf,"section %s not processed",name);
		  print_error(WARNING,"read_all",erbuf,NULL);
		}
	      ret=fscanf(inp_ops_fil,"%499[^\n\r]%*[\n\r]",buf_line);
	  }
    }
  return n;
} 

int read_hdr(char *fname, struct section *head)
{ // reads header, close stream on output
  int ret;
  char *sname="read_hdr";
  FILE* fil;
  fil=fopen(fname,"r");
  if(fil==NULL){
    print_error(FERROR,sname,"Can't open inp header",fname);
   
  }
  copy_ops=0;
  ret=read_all(fil,head,0); 
  if(ret<1)
    {
       print_error(FERROR,sname,"reading inp header",fname);
     
    }
    
  if(fclose(fil))
    print_error(WARNING,sname,"error closing inp header",fname);
    
  return ret;
}

int write_hdr(char *finp, char *fout, struct section *head)
{
  // copies header, closes input stream doesn't close output stream
  // header can not be pipe (it is read twice !)
  // if names NULL uses previously open files
 
  int ret;
  char *sname="write_hdr";
  if(finp)inp_ops_fil=fopen(finp,"r");
  if(inp_ops_fil==NULL){
     print_error(FERROR,sname,"Can't open inp header",finp);
   
  }
 
  if(fout)out_ops_fil=fopen(fout,"w");
  if(out_ops_fil==NULL){
     print_error(FERROR,sname,"Can't open out header",fout);
  }
  copy_ops=1;
  ret=read_all(inp_ops_fil,head,0);
  if(fclose(inp_ops_fil))
    print_error(WARNING,sname,"error closing inp header",finp);
  copy_ops=0;
  return ret;
}

char** cross_labels(char **chan_labels)
    {
	char **l,*b;
	int n,nn,i,j,k;
        if(chan_labels==NULL)
	   return NULL; 
	n=get_lab_no(chan_labels);
	nn=n*n;
	
	l=(char**)calloc(nn+1,sizeof(char*));
	if(!l)mem_error("cross_labels","l");
	set_lab_no(l,nn);
	k=1;
	for(i=0;i<n;i++)
	    {
		strncpy(buf_line,chan_labels[i+1],MAX_LINE/2-1);
		strcat(buf_line,"^");
		b=buf_line+strlen(buf_line);
		for(j=0;j<n;j++)
		    {
		    strncat(buf_line,chan_labels[j+1],MAX_LINE/2-1);
		    l[k++]=strdup(buf_line);
		    *b=0;
		    }
	    }
	return l;
    }
