/*************************************************
 *Copyright Piotr J. Franaszczuk 2006 *
 *ERL JHU
 **************************************************/
/*  procedures for data input and output for signal analysis */
/*  Piotr Franaszczuk JHU 2001,2002,2003 */ 
/*  corrected some read_next errors and piping */
/* added Start_chan to OUTPUT section */
/* no select if inp.dbuf==NULL   */

/* 12 nov 2002  chans_per_file==-1 if multiple files */
/* 16 dec 2002  output calibration=1 corrected */
/* 26 Jan 2003  nvars error in get_data corrected */
/* March 2003  modifications to get_data to speed up for no mtg */
/* July 2003   list file in memory to allow overlapped epochs */
/* October 2003 corrected time conversion for DST */
/* April 2004 added EDF read */
/* May 2004 added ID and  comments  to selection line and list line*/
/* Jan 2005 added next_file, next_selection, next_epoch,rewind_inp,rewind_all */
/* Jan 2005 added set_out_win, and option in init_io    54555454*/


#include "data_io.h"
const char dataio_ver[]="$Revision: 1.21 $";
extern char cfg_ver[];
extern char version[];
#ifdef WINDOWS
	#include <WinSock2.h>
#else
	extern long timezone;
#endif
char *progname;
char hst[100];
char *host=hst; 
char *ctl_name;
int merge_selections=0;  // default for get_epochs,allows overlapped

char *SPACE=" ";

/*-------------------------------------------------------------------*/	  
struct hfmt h_fmt={0}, v_fmt={0};
struct tag calib[]= {
    {"Chan_calib",   &h_fmt.chan_calib,calib_in, calib_out},
    {"Var_calib",&v_fmt.chan_calib,calib_in, calib_out}
};
struct section calib_sec={"CALIBRATION",calib,2,NULL};


/*-------------------------------------------------------------------*/	  

struct tag vars[]= {
    {"Var_labels",   &v_fmt.chan_labels,   labels_in,labels_out},
    {"AVG_definitions", &v_fmt.avg_defs, NULL,labels_out},// AVG definitions info  has to be vars[1]
    {"Var_units",&v_fmt.chan_units,labels_in,labels_out}
};
struct section vars_sec={"VARIABLES",vars,sizeof(vars)/sizeof(struct tag),&calib_sec};

/*-------------------------------------------------------------------*/	  

struct tag channels[]= {
    {"Chan_labels",   &h_fmt.chan_labels,   labels_in,labels_out},
    {"AVG_definitions", &h_fmt.avg_defs, NULL,labels_out},// AVG definitions info hast to be channels[1]
    {"Chan_units",&h_fmt.chan_units,labels_in,labels_out}
};
struct section channels_sec={"CHANNELS",channels,sizeof(channels)/sizeof(struct tag),&vars_sec};

/*-------------------------------------------------------------------*/	

struct tag  file_fmt[]= {
    {"Type", &h_fmt.title, str_in, str_out},  
    {"Samp_rate", &h_fmt.samp_rate, double_in, double_out},
    {"Bytes_per_win", &h_fmt.bytes_per_win, int_in,   int_out},
    {"Numb_chans",&h_fmt.nchans,int_in,   int_out},
    {"Byte_order",&h_fmt.byte_order,enum_in,  enum_out},
    {"File_format",   &h_fmt.file_fmt,  enum_in,  enum_out},
    {"Bytes_per_var", &h_fmt.bytes_per_var, int_in, int_out},
    {"Numb_vars", &h_fmt.no_of_vars,int_in, int_out},
    {"Data_offset",   &h_fmt.offset,int_in,   int_out},
    {"Chans_per_file",&h_fmt.chans_per_file, int_in,   int_out},
    {"Win_len",   &h_fmt.win_len,   double_in,  double_out}, 
    {"Win_shift", &h_fmt.win_shift, double_in,  double_out},
    {"Win_rate",  &h_fmt.win_rate,  double_in,  double_out},
    {"Name_template", &h_fmt.templ, str_in, str_out},
    {"Start_chan_no", &h_fmt.fil_ch_no, int_in, int_out},
    {"List_file", &h_fmt.lst_name,  str_in, str_out},
    {"Varmtg_file",   &h_fmt.vmtg_name, str_in, str_out}, 
    {"Montage_file",  &h_fmt.mtg_name,  str_in, str_out} 
};
struct section file_fmt_sec=
{"FILE_FMT",file_fmt,sizeof(file_fmt)/sizeof(struct tag),&channels_sec};

/*-------------------------------------------------------------------*/	 
struct pinfo p_info;
struct tag patient_info[]={
    {  "ID",&p_info.ID,   str_in,   str_out},
    {"First_Name",&p_info.fname,str_in,   str_out},
    {"Last_Name", &p_info.lname,str_in,   str_out},
    {"DOB",   &p_info.dob,  str_in,   str_out}
}; 
struct section patient_info_sec={
    "PATIENT_INFO",patient_info,sizeof(patient_info)/sizeof(struct tag),&file_fmt_sec};



/*-------------------------------------------------------------------*/	
struct inpdef inp_dflt={0},inp_def={0}, *in=&inp_def; 
struct tag data_inp[]={
    {"Type",   &inp_def.title,   str_in, NULL},
    {"Path",&inp_def.path,   str_in, NULL},
    {"Header_file", &inp_def.hdr_name,   str_in, NULL},
    {"Select_file", &inp_def.sel_name,   str_in, NULL},
    {"Montage_file",&inp_def.mtg_name,   str_in, NULL},
    {"Varmtg_file", &inp_def.vmtg_name,  str_in, NULL},
    {"Epoch",   &inp_def.epoch,  str_in, NULL},
    {"Win_len", &inp_def.win_len,double_in, NULL},
    {"Win_shift",   &inp_def.win_shift,  double_in, NULL},
    {"Calibrate",   &inp_def.calib,  int_in, NULL},  
    {"Numb_points", &inp_def.npoints,int_in, NULL},
    {"ID_select",   &inp_def.id_select,  int_in, NULL,},
    {"Shift_points",&inp_def.spoints,int_in, NULL},
    {"Put_zeros",   &inp_def.zeros,int_in,NULL}
};
struct section data_inp_sec={
    "INPUT", data_inp,sizeof(data_inp)/sizeof(struct tag),NULL};

/*-------------------------------------------------------------------*/	 
struct outdef out_dflt={0},out_def={0},*ou=&out_def;
struct tag data_out[]={
    {"Type", &out_def.title,str_in, NULL},
    {"Path",  &out_def.path, str_in, NULL},
    {"Header_file",   &out_def.hdr_name, str_in, NULL},
    {"List_file", &out_def.lst_name, str_in, NULL},
    {"Name_template", &out_def.templ,str_in, NULL},
    {"Byte_order",&out_def.byte_order, enum_in, NULL},
    {"File_format",   &out_def.file_fmt,   enum_in, NULL},
    {"Numb_vars", &out_def.no_of_vars, int_in, NULL},
    {"Start_chan",&out_def.start_chan, int_in,NULL},
    {"Numb_chans",&out_def.no_of_chans,int_in, NULL},
    {"All_chans", &out_def.all_chans,  int_in, NULL},
    {"Max_len",   &out_def.max_len,int_in, NULL},
    {"Input_no",  &out_def.inp_no, int_in, NULL},
    {"Win_len",   &out_def.win_len,  double_in, NULL},
    {"Win_shift", &out_def.win_shift,double_in, NULL},
    {"Calibrate", &out_def.calib,int_in, NULL},  
    {"Alarm", &out_def.alarm,int_in, NULL},
    {"Samp_rate", &out_def.samp_rate,  double_in, NULL}, 
    {"Start_chan_no", &out_def.fil_ch_no,int_in,   NULL},
    {"Chan_labels",   &out_def.chan_labels,labels_in, NULL},
    {"Chan_units",&out_def.chan_units, labels_in, NULL},
    {"Var_labels",&out_def.var_labels,labels_in, NULL},
    {"Var_units", &out_def.var_units, labels_in, NULL},
    {"Join",&out_def.join,int_in,NULL},
    {"Chans_per_file",&out_def.chans_per_file, int_in, NULL}
};
struct section data_out_sec={
    "OUTPUT", data_out,sizeof(data_out)/sizeof(struct tag),&data_inp_sec};
int inp_mode=mode_parallel; // default
int epoch_def=selections;
/*-------------------------------------------------------------------*/	 

int no_of_inp=-1,no_of_out=-1;   // new default is negative i.e. counting sections
struct tag inp_out_def[]={
    {"Numb_inputs", &no_of_inp, int_in, NULL},
    {"Mode",&inp_mode,  enum_in, NULL},
    {"Epochs",&epoch_def,enum_in,NULL},
    {"Merge_selections",&merge_selections,int_in,NULL},
    {"Numb_outputs",&no_of_out, int_in, NULL}
};
struct section inp_out_sec={
    "INPUT_OUTPUT", inp_out_def,sizeof(inp_out_def)/sizeof(struct tag),&data_out_sec};

/*-------------------------------------------------------------------*/	 

// this is only for parallel reading it can't be combined now with sequential reading
char **all_labels, **all_units;
int all_chans, parallel_read=1;  
double beg_first;// beg of first list needed by check_list 
char *first_title;
int cur_inp=0,inp_eof=0, inp_sel=0;


struct section *hdr_sec_list=&patient_info_sec;
struct section *io_sec_list=&inp_out_sec;
//struct section *last_sec =&calib_sec;

// all externs defined in cfg_io.c 
extern char *del_list[];  //max no of sections to delete
extern int del_no;
extern char *delimiters;   /* has to be global */
extern FILE *inp_ops_fil, *out_ops_fil;
extern char buf_line[MAX_LINE];
extern int copy_ops;/* default no copy in read_all */

int put_zeros(struct inpdef *inp,void *buf,int m,int n)
{ 
	int k=inp->end_pos-inp->pos;
	if(k>n)k=n;
	bzero(buf,m*k);
	return k;
}


int irint(double x)
{
    // replaces (int)rint
    // takes care of overflows in coversion from double to int
    if(x<-(double)INT_MAX)
	return -INT_MAX;
    else if(x>(double)INT_MAX)
	return INT_MAX;
    else
	return (int)rint(rint(10.*x)/10.); // to avoid problems close to 0.5
}

  

int prec(double x)
{
// counts number of trailing zero bits in mantissa of x
// it does not always tell about real precision !!! 
    union {double x; int n[2];}u;
    unsigned int b;
    int count;
    u.x=x;
#if(MY_BYTE_ORDER==little_endian)
    b=u.n[0];
    u.n[0]=u.n[1];
    u.n[1]=b;
    swab(&u,&u,4); 
#endif
    count=0;
    if((b=u.n[1]))
    {
	while((b&1))
	{
	    b>>=1;
	    count++;
	}
	return count;
    }
    else
	count=32;

    if((b=u.n[0]))
    {
	while(!(b&1))
	{
	    b>>=1;
	    count++;
	}
	return count;
    }
    else
	count=0; // both zero only for x=0
    return count;
}

	 	
struct inpdef* get_inp_str(struct inpdef* start,char* title)
{
    // returns inp stream with  title
    // if not found returns NULL
    // starts looking after stream start

    int i,is;
    if(start)
	is=start-in+1;
    else
	is=0;
    for(i=is;i<no_of_inp;i++)
    {
	if(strcmp(title,in[i].title))continue;
	return in+i;
    }
    return NULL;
}
struct outdef* get_out_str(struct outdef* start,char* title)
{
    // returns out stream with  title
    // if not found returns NULL
    // starts looking after stream start
   // or from begining if start=NULL

    int i,is;
    if(start)
	is=start-ou+1;
    else
	is=0;
    for(i=is;i<no_of_out;i++)
    {
	if(strcmp(title,ou[i].title))continue;
	return ou+i;
    }
    return NULL;
}


void write_out_hdr(struct outdef *out)
{
    // write out header based on out  values

    char name[MAX_FILE_NAME],name1[MAX_FILE_NAME],*sname="write_out_hdr";
    struct inpdef *inp;
    int i;

    h_fmt.samp_rate=out->samp_rate;  		 
    h_fmt.win_len=out->win_len;
    h_fmt.win_shift=out->win_shift;
    if(out->rate_flag)
	h_fmt.win_rate=out->win_rate;
    else
	h_fmt.win_rate=-1.;
	 
    h_fmt.file_fmt=out->file_fmt;
    h_fmt.bytes_per_win=out->bytes_per_win;
    h_fmt.byte_order=out->byte_order;
    h_fmt.chans_per_file=out->chans_per_file;
    h_fmt.fil_ch_no=out->fil_ch_no;
    h_fmt.nchans=out->no_of_chans;
    if(h_fmt.title)free(h_fmt.title);
    h_fmt.title=out->title;
    if(out->var_len)
	h_fmt.no_of_vars=-out->no_of_vars;
    else
	h_fmt.no_of_vars=out->no_of_vars;

    h_fmt.bytes_per_var=out->bytes_per_var;

    if(h_fmt.no_of_vars>0)
	h_fmt.bytes_per_win=h_fmt.no_of_vars* h_fmt.bytes_per_var* h_fmt.nchans;
    if(h_fmt.lst_name)free(h_fmt.lst_name);
    h_fmt.lst_name=out->lst_name;
    h_fmt.offset=0;
    if(out->ext)
	sprintf(name,"%s#%0*d.%s",out->templ,out->templ_dig,out->file_no,out->ext);
    else
	sprintf(name,"%s#%0*d",out->templ,out->templ_dig,out->file_no);
    if(h_fmt.templ)free(h_fmt.templ);
    h_fmt.templ=strdup(name);
    if(h_fmt.chan_labels)
    {
	for(i=0;i<get_lab_no(h_fmt.chan_labels);i++)
	    free(h_fmt.chan_labels[i+1]);
	free(h_fmt.chan_labels);
 
    }
    h_fmt.chan_labels=out->chan_labels;
    if(h_fmt.chan_units)
    {
	for(i=0;i<get_lab_no(h_fmt.chan_units);i++)
	    free(h_fmt.chan_units[i+1]);
	free(h_fmt.chan_units);
    }
    h_fmt.chan_units=out->chan_units;

    if(v_fmt.chan_labels)
    {
	for(i=0;i<get_lab_no(v_fmt.chan_labels);i++)
	    free(v_fmt.chan_labels[i+1]);
	free(v_fmt.chan_labels);
    }
    v_fmt.chan_labels=out->var_labels;
    if(v_fmt.chan_units)
    {
	for(i=0;i<get_lab_no(v_fmt.chan_units);i++)
	    free(v_fmt.chan_units[i+1]);
	free(v_fmt.chan_units);
    }
    v_fmt.chan_units=out->var_units;

 

    if(out->inp_no<0 || out->inp_no>no_of_inp-1)
    {
	print_error(WARNING,sname,"wrong input_no  or  INPUT section not defined, set to 0",NULL);
	out->inp_no=0; 
    }
 
    if(*out->hdr_name=='/'||!strncmp(out->hdr_name,"./",2))
	strcpy(name,out->hdr_name);
    else
    {
	strcpy(name,out->path);
	strcat(name,out->hdr_name);
    } 
    del_list[del_no++]="CALIBRATION";
    if(no_of_inp>0)
    { 
	inp=in+out->inp_no;
	if(*inp->hdr_name=='/'||!strncmp(inp->hdr_name,"./",2))
	    strcpy(name1,inp->hdr_name);
	else
	{
	    strcpy(name1,inp->path);
	    strcat(name1,inp->hdr_name);
	}
  
	if(inp->mtg.navg<3)
	    channels[1].fout=NULL;
 
	if(inp->vmtg.navg<32)
	    vars[1].fout=NULL;

	if(inp->mtg_name)
	{
	    if(h_fmt.mtg_name)free(h_fmt.mtg_name);
	    h_fmt.mtg_name=inp->mtg_name;
	} 
	write_hdr(name1,name,hdr_sec_list);
	out->hdr=out_ops_fil;  //save it for use outside
	if(out->calib)
	{
	    // keep calibration from mtg
	    del_list[--del_no]=NULL;
	    if(inp->mtg.chan_calib)
	    {
		if(!(h_fmt.chan_calib=(double*)calloc(out->no_of_chans+1,sizeof(double))))
		    mem_error(sname,"chan_calib");
		h_fmt.chan_calib[0]=out->no_of_chans;
		for(i=1;i<=out->no_of_chans;i++)
		    h_fmt.chan_calib[i]=inp->mtg.chan_calib[i+out->start_chan];
	    }
	    if(inp->vmtg.chan_calib)
	    {
		// var calibration
		if(!(v_fmt.chan_calib=(double*)calloc(out->no_of_vars+1,sizeof(double))))
		    mem_error(sname,"var_calib");
		v_fmt.chan_calib[0]=out->no_of_vars;
		for(i=1;i<=out->no_of_vars;i++)
		    v_fmt.chan_calib[i]=inp->vmtg.chan_calib[i+out->start_chan];
	    }

	    write_all(out->hdr,&calib_sec);
	}
    }
    else
    {   
	if(!(out->hdr=fopen(name,"w")))
	    print_error(FERROR,sname,"Can't open out header",name);
	out_ops_fil=out->hdr;
	write_all(out->hdr,hdr_sec_list);
	if(out->calib)
	    // calibration should be already in h_fmt
	    write_all(out->hdr,&calib_sec); 
	 
    }

    channels[1].fout=labels_out;
    vars[1].fout=labels_out;


    // reset header
    bzero(&h_fmt,sizeof(h_fmt));
    bzero(&v_fmt,sizeof(v_fmt));
    fflush(out->hdr);
}
	
	
int read_mtg(char* fname, struct mont *mtg, struct hfmt *fmt )
{ /* reads montage and puts it in  structure mtg */
  /* it has to be called after reading header to h_fmt */

    FILE* fil;
    int i,j,n,k,av_defs;
    char *s,*a,*sname="read_mtg";
    int *tmp;
    char err_buf[500];

    /*
      if(fmt->chan_labels==NULL)
      {
      if(fmt==&h_fmt)a="C";
      else a="V";
      fmt->chan_labels=malloc((fmt->nchans+1)*sizeof(char*));
      set_lab_no(fmt->chan_labels,fmt->nchans);
      for(i=1;i<=fmt->nchans;i++)
      {
      char buf[MAX_LINE];
      sprintf(buf,"%s%d",a,i);
      fmt->chan_labels[i]=strdup(buf);
      }
	
      }
    */
    if(fname==NULL)
    { 

	/* copy from header all channels info */
	/* ignore empty labels marked by . */ 
	/* for variables it is not needed */  


	mtg->type=NO_MTG;  // no change at all
	mtg->navg=0;
	mtg->lst_avg=NULL;
	mtg->ref_no=mtg->chan_no=mtg->ch_no=NULL;
	if(fmt->chan_labels && fmt->nchans<1000)
	{
	    for(i=0,k=0;i<fmt->nchans;i++)
		if(strlen(fmt->chan_labels[i+1])>0 && fmt->chan_labels[i+1][0]!='.')k++;
	    mtg->nchans=mtg->mchans=k;
	    if(!(mtg->chan_labels=(char**)malloc((mtg->nchans+1)*sizeof(char*))))
		mem_error(sname,"mtg.chan_labels");
	    set_lab_no(mtg->chan_labels,mtg->mchans);
	    if(fmt==&h_fmt)
	    {
		if(!(mtg->chan_no=(int*)malloc(mtg->nchans*sizeof(int))))
		    mem_error(sname,"mtg.chan_no");
 
		if(!(mtg->ch_no=(int*)malloc(mtg->mchans*sizeof(int))))
		    mem_error(sname,"mtg.ch_no");
		for(i=0,k=0;i<fmt->nchans;i++)
		    if(strlen(fmt->chan_labels[i+1])>0 && fmt->chan_labels[i+1][0]!='.')
		    { 
			mtg->chan_labels[k+1]=strdup(fmt->chan_labels[i+1]);
			mtg->ch_no[k]=i;
			mtg->chan_no[k]=k;k++;
		    } 
	    }
	}
	else
	{
	    mtg->chan_labels=NULL;
	    mtg->nchans=mtg->mchans=fmt->nchans;
        }
            if(fmt==&h_fmt)
            {
	       if(!(mtg->chan_no=(int*)malloc(mtg->nchans*sizeof(int))))
		    mem_error(sname,"mtg.chan_no");
 
		if(!(mtg->ch_no=(int*)malloc(mtg->mchans*sizeof(int))))
		    mem_error(sname,"mtg.ch_no");
		for(i=0;i<fmt->nchans;i++)
               		mtg->ch_no[i]=mtg->chan_no[i]=i;
             }
	if(mtg->nchans!=fmt->nchans)
	    mtg->type=MTG_SEL;
 
    }
    else
    {
	if(!fmt->chan_labels)
	    print_error(FERROR,sname,"No header labels for montage",fname);
	mtg->type=MTG_SEL;  // no change in ref but may change in channel no
	fil=fopen(fname,"r");
	if(fil==NULL)
	    print_error(FERROR,sname,"Can't open for read",fname);
	fscanf(fil,"%d\n",&mtg->nchans);
	n=mtg->nchans;
	if(n<1)
	    print_error(WARNING,sname,"No_of_channels<1",fname);

	if(!(mtg->chan_labels=(char**)malloc((mtg->nchans+1)*sizeof(char*))))
	    mem_error(sname,"mtg.chan_labels");
 
	if(!(mtg->chan_no=(int*)malloc(n*sizeof(int))))
	    mem_error(sname,"mtg.chan_no");
	if(!(mtg->ref_no=(int*)malloc(n*sizeof(int))))
	    mem_error(sname,"mtg.ref_no");
 
	if(!(tmp=(int*)calloc(fmt->nchans+1,sizeof(int))))
	    mem_error(sname,"tmp_chan");
	mtg->navg=0;  
 
	set_lab_no(mtg->chan_labels,n);
	for(i=0;i<n;i++)
	{  //read mtg labels and defs  
	    if(fscanf(fil,"%99[^\n]\n",buf_line)!=1)
	   
	    {
		sprintf(err_buf,"error in line %d:%s\n",i+2,buf_line);
		print_error(FERROR,sname,err_buf,fname);
	    }
	    if(buf_line[0]!='$')
	    { char *c;
	    s=buf_line;
	    if(!(mtg->chan_labels[i+1]=strdup(s)))mem_error(sname,"mtg.chan_labels[i]"); 
	    c=strchr(mtg->chan_labels[i+1],'-');
	    if(c)*c='~';
	    }  
	    else
	    {
		s=strtok(buf_line+1,"=");
		if(!(mtg->chan_labels[i+1]=strdup(s)))mem_error(sname,"mtg.chan_labels[i]");
		s=NULL;
	    }
	    a=strtok(s,"-");
	    if(a==NULL) 
	    {
	  
		print_error(FERROR,sname,"Empty chan label",fname);
	    }

	    for(j=0;j<fmt->nchans;j++)
//		if(!strncmp(a,fmt->chan_labels[j+1],strlen(a)))break;
	    if(!strcmp(a,fmt->chan_labels[j+1]))break;
	    if(j==fmt->nchans && strncmp(a,"AVG",3))
	    {
		sprintf(err_buf,"chan label(1) %s not recognized",a);
		print_error(FERROR,sname,err_buf,fname);
	    }
	    else if(j<fmt->nchans)
	    {
		mtg->chan_no[i]=j;
		tmp[j+1]=1;
	    }
	    else
	    {	
		// average
		k=atoi(a+3);
		if(k+1>mtg->navg)
		    mtg->navg=k+1;
		mtg->chan_no[i]=j+k; 
		// AVG0 is nchans and so on   
	    }
 
	    a=strtok(NULL,"-");
	    if(a==NULL)
	    {
		mtg->ref_no[i]=-1;
		tmp[0]=1;
	    }
	    else
	    {
		mtg->type=MTG_REF; // simple subtraction in at least one chan
		for(j=0;j<fmt->nchans;j++)
		    if(!strncmp(a,fmt->chan_labels[j+1],strlen(a)))break;
		if(j==fmt->nchans && strncmp(a,"AVG",3))
		{
		    sprintf(err_buf,"chan label(2) %s not recognized",a);
		    print_error(FERROR,sname,err_buf,fname);		  
		}
		else if(j<fmt->nchans)
		{
		    mtg->ref_no[i]=j;
		    tmp[j+1]=1;
		}
		else
		{	
		    // average
		    k=atoi(a+3);
		    if(k+1>mtg->navg)
			mtg->navg=k+1;
		    mtg->ref_no[i]=j+k; 
		    // AVG0 is nchans and so on
		}
	    }
	}// end for loop 
  
	if(mtg->navg)
	{ // definitions of averages
	  // AVG0 predefined as average of all in input
	  // AVG1 predefined as average of all present in mtg
	  // rest has to be defined by user in mtg_file
	    mtg->type=MTG_AVR ;  // at least one channel has average ref	  
	    if(!(mtg->lst_avg=(int**)calloc(mtg->navg,sizeof(int*))))
		mem_error(sname,"mtg.lst_avg");
	    if(!(mtg->avg_no=(int*)calloc(mtg->navg,sizeof(int))))
		mem_error(sname,"mtg.lst_avg");
	    // find which averages are needed
	    for(i=0;i<mtg->nchans;i++)
	    {
		if((k=mtg->chan_no[i]-fmt->nchans)>=0)mtg->avg_no[k]++;  //temporary >0
		if((k=mtg->ref_no[i]-fmt->nchans)>=0)mtg->avg_no[k]++;
	    }

	    if(mtg->avg_no[1]>0)
	    {
		// count how many different channels for AVG1
		for(i=1,j=0;i<=fmt->nchans;j+=tmp[i++]);
  
		mtg->mchans=mtg->avg_no[1]=j;
		if(!(mtg->lst_avg[1]=malloc(mtg->mchans*sizeof(int))))
		    mem_error(sname,"mtg.lst_avg[1]");
		for(i=1,j=0;i<=fmt->nchans;i++)
		    if(tmp[i])mtg->lst_avg[1][j++]=i-1; //original channel no
	    }
	    // it has to be in this order (first 1 and 0 )
	    if(mtg->avg_no[0]>0)
	    {
		mtg->mchans=fmt->nchans; //for AVG0
		for(i=0;i<=fmt->nchans;tmp[i]=i++);

		mtg->avg_no[0]=mtg->mchans;
		mtg->lst_avg[0]=NULL;  //no nead for list
	    }
	    av_defs=0;
	    for(i=2;i<mtg->navg;i++)
	    {
		//	  if(mtg->avg_no[i]>0)
		{
		    av_defs+=mtg->avg_no[i];
		    if((fscanf(fil,"%99[^\n]\n",buf_line)!=1) || 
		       (sscanf(buf_line,"AVG%d=%d\n",&j,&n))!=2 || n<=0)
		    {
			sprintf(err_buf,"expected def of AVG%d in line %s",i,buf_line);
			print_error(FERROR,sname,err_buf,fname);
		    }
		    if(!(mtg->lst_avg[i]=malloc(n*sizeof(int))))
			mem_error(sname,"mtg.lst_avg[i]");
		    if(mtg->avg_no[i]>0)
			mtg->avg_no[i]=n;
		  
		    for(j=0;j<n;j++)
		    { int jj;
		    if(fscanf(fil,"%99[^\n]\n",buf_line)!=1)
		    {
			sprintf(err_buf,"in def of AVG%d label %d:%s",i,j,buf_line);
			print_error(FERROR,sname,err_buf,fname);
		    }  
		    if(!mtg->avg_no[i])continue;// skipping definition if never used
		    a=buf_line;
		    for(jj=0;jj<fmt->nchans;jj++)
			if(!strncmp(a,fmt->chan_labels[jj+1],strlen(a)))break;
		  
		    if(jj==fmt->nchans && strncmp(a,"AVG",3))
		    {
			sprintf(err_buf,"chan label %s in def of AVG%d not recognized",a,i);
			print_error(FERROR,sname,err_buf,fname);
		    }
		    if(jj<fmt->nchans)
		    {
			tmp[jj+1]=1;
			mtg->lst_avg[i][j]=jj;
		    }
		    else
		    {	
			// average
			k=atoi(a+3);
			if(k>=i)
			{
			    sprintf(err_buf,"definition of AVG%d can't include %s",i,a);
			    print_error(FERROR,sname,err_buf,fname);
			}
			else if(mtg->avg_no[k]==0)
			{
			    sprintf(err_buf,"AVG%d in def of AVG%d not defined",k,i);
			    print_error(FERROR,sname,err_buf,fname);
			}	
			mtg->lst_avg[i][j]=jj+k; 
		    }
		    }
		}
	  
	    } //end of for(i
	    if(av_defs)
	    {
		int jj,kk;
		av_defs+=mtg->navg-2;
		mtg->avg_defs=(char**)calloc(av_defs+1,sizeof(char*));
		set_lab_no(mtg->avg_defs,av_defs);
		for(i=2,jj=1;i<mtg->navg;i++)
		    if((n=mtg->avg_no[i])>0)
		    {
			sprintf(buf_line,"AVG%d=%d",i,n);
			mtg->avg_defs[jj++]=strdup(buf_line);
			for(j=0;j<n;j++)
			{
			    kk=mtg->lst_avg[i][j];
			    if(kk<fmt->nchans)
				mtg->avg_defs[jj++]=strdup(fmt->chan_labels[kk+1]);
			    else
			    {
				sprintf(buf_line,"AVG%d",kk-fmt->nchans);
				mtg->avg_defs[jj++]=strdup(buf_line);
			    }
			}	  
		    }
	    }
			  
 
	}// end of definitions of averages
	fclose(fil);


  	 
	if(mtg->navg!=1)
	{
	    // count how many different channels and assign numbers
	    j=0;
	    tmp[0]=j++;  //  dbuf[0]=0;
	    for(i=1;i<=fmt->nchans;i++)
		if(tmp[i])tmp[i]=j++;   
	    mtg->mchans=j-1;
	}
 
	if(!(mtg->ch_no=(int*)malloc(mtg->mchans*sizeof(int))))
	    mem_error(sname,"mtg.ch_no");
	if(mtg->navg==1)
	    for(i=0;i<fmt->nchans;i++)
		mtg->ch_no[i]=i;
	else
	    for(i=0,j=0;i<fmt->nchans;i++)
		if(tmp[i+1])mtg->ch_no[j++]=i;  //original chan_no
 
	for(i=0;i<mtg->nchans;i++)
	{
	    if(mtg->chan_no[i]<fmt->nchans)
		mtg->chan_no[i]=tmp[mtg->chan_no[i]+1];
	    else
		mtg->chan_no[i]+=mtg->mchans-fmt->nchans+1;
	    if(mtg->ref_no[i]<fmt->nchans)
		mtg->ref_no[i]=tmp[mtg->ref_no[i]+1];
	    else
		mtg->ref_no[i]+=mtg->mchans-fmt->nchans+1;
	}
   
	for(i=1;i<mtg->navg;i++)
	    for(j=0;j<mtg->avg_no[i];j++)
		if(mtg->lst_avg[i][j]<fmt->nchans)
		    mtg->lst_avg[i][j]=tmp[mtg->lst_avg[i][j]+1];
		else
		    mtg->lst_avg[i][j]+=mtg->mchans-fmt->nchans+1;
	free(tmp);
	if(mtg->type==MTG_SEL)
	{
	    free(mtg->ref_no);
	    mtg->ref_no=NULL;
	    for(i=0;i<mtg->mchans;i++)
		mtg->chan_no[i]--;
	}
	  
    } // end of NO_MTG if
 
    if(fmt->chan_calib && mtg->ch_no)
    {
	if(!(mtg->chan_calib=(double*)malloc((mtg->mchans+1)*sizeof(double))))
	    mem_error(sname,"mtg.chan_calib");
	mtg->chan_calib[0]=mtg->mchans;
	for(i=0;i<mtg->mchans;i++)
	    mtg->chan_calib[i+1]=fmt->chan_calib[mtg->ch_no[i]+1];
    }
    if(fmt->chan_units && mtg->ch_no)
    {
	if(!(mtg->chan_units=(char**)malloc((mtg->mchans+1)*sizeof(char*))))
	    mem_error(sname,"mtg.chan_units");
	set_lab_no(mtg->chan_units,mtg->mchans);
	for(i=0;i<mtg->mchans;i++)
	    if(!(mtg->chan_units[i+1]=strdup(fmt->chan_units[mtg->ch_no[i]+1])))
		mem_error(sname,"mtg.chan_unit[i]");  	   
    }

    return 0;
} 


double get_time(char *time, double tod)
{
    static char *month[12]={"Jan","Feb","Mar","Apr","May","Jun",
			    "Jul","Aug","Sep","Oct","Nov","Dec"};
    char mon[5],*sname="get_time";
    int m;
    time_t t;
    struct tm tim={0}; 

  
    sscanf(time,"%*s %d %3s %d",&tim.tm_mday,mon,&tim.tm_year);
    for(m=0;m<12;m++)
	if(!strncmp(mon,month[m],3))break;

    if(m<12)
	tim.tm_mon=m;
    else
    {
	sprintf(buf_line,"Wrong month %s in time info",mon);
	print_error(ERROR,sname,buf_line,NULL);
    }
    tim.tm_year-=1900;
    if((t= mktime(&tim))<0)
    {
 
	print_error(ERROR,sname," mktime failed time conversion",NULL);
    }
    t-=timezone;
    return (double)t+tod;
} 
  
 

int read_list_line(struct inpdef *inp)
{
    // reads one line from list_file
    // return 0 if ok or -1 if eof
    char time[30];
    double len,tod,last;
    char text[200],*sname="read_list_line";
    int id;
    do{
	id=0;
	*text=0;
	if( !fgets(buf_line,MAX_LINE,inp->lst) ||
	    sscanf(buf_line,"%s\t%29[^\t]\t%lf\t%lf\t%d\t%[^\n]\n",inp->name,time,&len,&tod,&id,text)<4)
	{
	    if(feof(inp->lst))return -1;
	    print_error(FERROR,sname,"reading from  list_file",inp->lst_name);
	}
    }while((buf_line[0]=='#' && !inp->zeros) || len<=0 || (inp->id_select>0 && id!=inp->id_select));
    last=inp->f_start;
    tod=rint(tod*inp->freq)/inp->freq;   //rounding to sampling rate
    len=rint(len*inp->freq)/inp->freq;   //rounding to sampling rate
    inp->f_start=get_time(time,tod);
    
    if(!inp->sel_comment)
    {
	inp->sel_id=id;
        if(strlen(text)<1)
         inp->sel_text=NULL;
	else
	 inp->sel_text=strdup(text);
    }

    if(inp->f_start<inp->f_end)
	print_error(WARNING,sname,"overlapping intervals in list_file",inp->lst_name);

    inp->f_end=inp->f_start+len;
    inp->cur+=irint((last-inp->f_start)*inp->freq);
    inp->end_pos=irint(len*inp->freq);
    inp->e_start_pos=irint((inp->e_start-inp->f_start)*inp->freq);
    if(inp->e_end<1e100)
      inp->e_end_pos=inp->e_start_pos+irint((inp->e_end-inp->e_start)*inp->freq);
    else
      inp->e_end_pos=INT_MAX;
  
    return 0;
}
  
int read_select_line(struct inpdef *inp)
{
    // reads one line from select file
    // return 0 if ok or -1 if eof
    // sets e_start and e_end

    char time[30],text[200],*sname="read_select_line";
    double len,tod;
    int id;
   
    if(merge_selections)
	{ int k=inp->cur_epoch;
	  if(k>=inp->list[0])return -1;
	  inp->sel_comment=0;
          inp->sel_out_id=inp->sel_id; // remember last for output
	  inp->sel_id=inp->id[k+1];
          inp->sel_out_text=inp->sel_text; // remember last for output
          if(inp->text[k+1])
	        inp->sel_text=strdup(inp->text[k+1]);
	  else
	    inp->sel_text=NULL;
	  inp->e_start=inp->list[2*k+1];
          inp->e_end=inp->list[2*k+2];
      }
   else if(!inp->sel)
	print_error(ERROR,sname,"NULL sel file handle",NULL);
   else
    {
	  
    do{
	id=0;
	*text=0;
	if(fscanf(inp->sel,"%29[^\t]\t%lf\t%lf\t%d\t%[^\n]\n",time,&len,&tod,&id,text)<3)
	{
	    if(feof(inp->sel))return -1;
	    print_error(FERROR,sname,"reading from  select_file",inp->lst_name);
	}
    }while(len<=0 || (inp->id_select>0 && id!=inp->id_select)|| (time[0]=='#'));
    inp->sel_out_id=inp->sel_id; 
    inp->sel_out_text=inp->sel_text;
    if(id)
    {
	inp->sel_id=id;
	if(text)
		inp->sel_text=strdup(text);
	else
	    inp->sel_text=NULL;
	inp->sel_comment=1;
    }
    else
	inp->sel_comment=0;
    tod=rint(tod*inp->freq)/inp->freq;   //rounding to sampling rate
    len=rint(len*inp->freq)/inp->freq;   //rounding to sampling rate
    inp->e_start=get_time(time,tod);
//    xx=inp->e_start-inp->f_start;
    inp->e_end=inp->e_start+len;
 //   printf("e_start=%lf,f_start=%lf,freq=%lf\n",inp->e_start,inp->f_start,inp->freq);
    }
 //   inp->e_start_pos=irint(xx*inp->freq);
   inp->e_start_pos=irint((inp->e_start-inp->f_start)*inp->freq);
  
   if(inp->e_end<1e100)
      inp->e_end_pos=inp->e_start_pos+irint((inp->e_end-inp->e_start)*inp->freq);
    else
      inp->e_end_pos=INT_MAX;  
 
    inp->new_epoch=1;
    //if(len<=0)
    //print_error(WARNING,sname,"last<0",NULL);
    return 0;
}
  

int write_list_line(struct outdef *out)
{
    // writes one line to list_file
    // uses values from global vars out_start and out_len
    // if exiast puts out->sel_id and out->sel_text

    char time[30],*sname="write_list_line";
    double tod;
    time_t t;
    struct tm *tim;
    int id=out->sel_id;
    char* text=out->sel_text;
    out->start=rint(out->start/out->win_shift)*out->win_shift;   //rounding to sampling rate
    out->len=rint(out->len/out->win_shift)*out->win_shift;   //rounding to sampling rate
    t=(time_t)(out->start);
    tim=gmtime(&t);

    strftime(time,30,"%a %d %b %Y %H:%M",tim); 
    tim->tm_hour=tim->tm_min=tim->tm_sec=0;
    tod=out->start-mktime(tim)+timezone;
    
    if(fprintf(out->lst,"%s\t%s\t%.4lf\t%.4lf",out->name,time,out->len,tod)<1)
	goto error;
    if(id>0 && fprintf(out->lst,"\t%d",id)<1)
	goto error;
    if(text && *text && fprintf(out->lst,"\t%s",text)<1)
	goto error;
    if(fprintf(out->lst,"\n")<1)
	goto error;

    fflush(out->lst);
  
    return 0;
 error:print_error(FERROR,sname,"writing to list_file",out->lst_name);
  
}

int write_select_line(FILE* sel,double dt,double t_cur,double len,int no,char*type)
{
    // writes one line to select file sel
    // t_cur is begining len length
 
 
    char time[30],*sname="write_select_line";
    double tod;
    time_t t;
    struct tm *tim;
    //  fprintf(stderr,"%s:in write_list_line\n",progname);
    //  t_cur=inp->f_cur+offset;
    if(len<0)
    {
	t_cur+=len;
	len=-len;
    }
    t_cur=rint(t_cur/dt)*dt;   //rounding to sampling rate
    len=rint(len/dt)*dt;   //rounding to sampling rate
    t=(time_t)(t_cur);
    tim=gmtime(&t);

    strftime(time,30,"%a %e %b %Y %H:%M",tim); 
    tim->tm_hour=tim->tm_min=tim->tm_sec=0;
    tod=t_cur-mktime(tim)+timezone;

    if(no)
    {
	if(fprintf(sel,"%s\t%.4lf\t%.4lf\t%d\t%s\n",time,len,tod,no,type)<1)
	    print_error(ERROR,sname,"writing to select_file",NULL);  
    }
    else
    {

	if(fprintf(sel,"%s\t%.4lf\t%.4lf\n",time,len,tod)<1)
	    print_error(ERROR,sname,"writing to select_file",NULL); 
    } 

    fflush(sel);
    //  fprintf(stderr,"%s:out write_list_line\n",progname);
    return 0;
}

void init_new_out(struct outdef *out,FILE* fil)
{ // generates new out file and writes to list file
  // fil is pointer to last file
  // return pointer to new name
  // writes to list
  // at start it should have out->start set or if <0. it is taken from inp->f_cur
  // on return out->start is set to inp->f_cur and out->count=0
  // out->new is ==1 on return if started new file
  // sets out->name
 
    
    char name[MAX_FILE_NAME],*sname="init_new_out";

    if(!parallel_read)
	out->inp_no=cur_inp;
    if(fil)
    {
	fclose(fil);
	fil=NULL; 
	if(out->rate_flag)
	    out->len=(out->count-1)/out->win_rate+out->win_len;
	else
	    out->len=(out->count-1)*out->win_shift+out->win_len;

	if(out->lst && !(out->pipe))
	    write_list_line(out);  

    }
    else
	if(out->inp_no>-1 && out->start<0.)  // needed for first
	{
	    out->start=in[out->inp_no].f_cur+out->offset; 
	    out->sel_id=in[out->inp_no].sel_out_id;
            if(out->sel_text)free(out->sel_text);
	    if(in[out->inp_no].sel_out_text)
		out->sel_text=strdup(in[out->inp_no].sel_out_text);
	    else
	       out->sel_text=NULL;
        }

    if(!out->lst)
    { // first write to lst; hdr should be closed by then for pipes to work
	if(out->hdr && fclose(out->hdr))
	{
	    print_error(WARNING,sname,"error closing hdr file",out->hdr_name);
	    out->hdr=NULL;
	}

	if(*out->lst_name=='/'||!strncmp(out->lst_name,"./",2))
	    strcpy(name,out->lst_name);
	else
	{
	    strcpy(name,out->path);
	    strcat(name,out->lst_name);
	}
	out->lst=fopen(name,"w");   
	if(!out->lst)
	{ 
	    print_error(FERROR,sname,"Can't open output list file",name);
	}
    }

    if(out->ext)
	sprintf(out->name,"%s%0*d.%s",out->templ,out->templ_dig,out->file_no++,out->ext);
    else
	sprintf(out->name,"%s%0*d",out->templ,out->templ_dig,out->file_no++);


 
    if(out->pipe) // for pipes write at the begining 
    {  
	if(out->start<0. && out->inp_no>-1)  // needed for first
	{
	    out->start=in[out->inp_no].f_cur+out->offset; 
	    out->sel_id=in[out->inp_no].sel_id;
	    if(out->sel_text)free(out->sel_text);
		 if(in[out->inp_no].sel_text)
		out->sel_text=strdup(in[out->inp_no].sel_text);
	    else
	       out->sel_text=NULL;
        }
	    out->len=1.e100;  // infinity
 
	write_list_line(out);
 
    }
    out->new_out=1;
    if(out->pipe)
	out->start=-1.;  //dflt for pipes
    else if (out->inp_no>-1)
    {
	out->start=in[out->inp_no].f_cur+out->offset;  // dflt start of next one
	out->sel_id=in[out->inp_no].sel_out_id;
	if(out->sel_text)free(out->sel_text);
	if(in[out->inp_no].sel_out_text && strlen(in[out->inp_no].sel_out_text)>0) out->sel_text=strdup(in[out->inp_no].sel_out_text);
        else out->sel_text=NULL;
    }
    else
	out->start+=out->len;

}
void new_out_file(struct outdef *out)
{ 
    // opens new file, only for single file
    // if (out->join) returns without opening 
    char name[MAX_FILE_NAME],*sname="new_out_file";
  
    if(out->join)return;
    init_new_out(out,out->fil);
    strcpy(name,out->path);
    strcat(name,out->name);
#ifndef _WIN32
    if(out->pipe)
    {
	unlink(name);   // delete pipe first
	if(mkfifo(name,S_IRUSR|S_IWUSR))
	{
	    print_error(FERROR,sname,"Can't create out fifo pipe",name);
	}
    }
#endif
    if(!(out->fil=fopen(name,"wb")))
    {
	print_error(FERROR,sname,"Can't open output file",name); 
    }
    out->count=out->bytes=0;
}



void new_out_file_chan(struct outdef *out,int chan)
{ 
    // if chan==0 call init_out
    // opens file for chan
    char name[MAX_FILE_NAME],*sname="new_out_file_chan"; 
    if(out->join)return;
    if(chan==0)
	init_new_out(out,out->files[0]);
    else
	if(out->files[chan]) fclose(out->files[chan]);

    sprintf(name,"%s%s.%3.3d",out->path,out->name,chan+out->fil_ch_no);
#ifndef _WIN32
    if(out->pipe)
    {
	unlink(name);
	if(mkfifo(name,S_IRUSR|S_IWUSR))
	{
	    print_error(FERROR,sname,"Can't create out fifo pipe",name); 
	}
    }
#endif
    if(!(out->files[chan]=fopen(name,"wb")))
    {
	print_error(FERROR,sname,"Can't open output file",name); 
    } 
    out->mbytes[chan]=0;
    if(chan==0)
	out->count=0;
    if(chan==(out->no_of_chans-1))
	out->bytes=0;
}
void close_inp_file (FILE* fil,char *name)
{
    if( fil!=NULL)
    {
	fclose(fil);
	fil=NULL;
	if(name[0])
	{	  
	    if(unlink(name))  // delete pipe
		print_error(WARNING,"close_inp_file", "Can't unlink pipe",name);
	}  
    }
}

  
  
int new_inp_file(struct inpdef *inp, int flag)
{ // return -1 for eof 0 if OK
  //  reads new line from list if flag (for pipes it should always be 1)
  // resets pos
    char name[MAX_FILE_NAME],*sname="new_inp_file"; 
    int i;
    if(inp->files)
	for(i=0;i<inp->nchans;i++)
	{ 
	    if(inp->files[i])
	    {
		fclose(inp->files[i]);
		inp->files[i]=NULL;
		if(inp->pipe)
		{
		    sprintf(name,"%s%s.%3.3d",inp->path,inp->name,inp->fil_ch_no+inp->ch_no[i]);
		    if(unlink(name))  // delete pipe
			print_error(WARNING,sname, "Can't unlink pipe",name);
		}
	    }
	}
    else if(inp->fil)
    {
	fclose(inp->fil);
	inp->fil=NULL;
	if(inp->pipe)
	{
	    strcpy(name,inp->path);
	    strcat(name,inp->name);
	    if(unlink(name))  // delete pipe
		print_error(WARNING,sname, "Can't unlink pipe",name);
	}
    }
 
    if(flag && read_list_line(inp))return -1;
    inp->zeros=(inp->name[0]=='#')?1:0;
    if(inp->files)
    {
	if(!inp->zeros)
	{
	    for(i=0;i<inp->nchans;i++)
	    {
	  
		sprintf(name,"%s%s.%3.3d",inp->path,inp->name,inp->fil_ch_no+inp->ch_no[i]);

		if(!(inp->files[i]=fopen(name,"rb")))
		{
#ifndef _WIN32
		    if(inp->pipe)
		    {
			int count=5; // wait 5 secs for pipe to open
			while(errno==ENOENT && count--)
			{
			    sleep(1); 
			    if((inp->files[i]=fopen(name,"rb")))goto open1; 
			}
		    }
#endif
	 
		    print_error(FERROR,sname,"Can't open input file",name);
	 
		}
	    open1:;
	    }
	    if(inp->offset)
	    {inp->pos=-1;skip_fix_many(inp,0);}  
	}
    }
    else
    {
	if(!inp->zeros)
	{
	strcpy(name,inp->path);
	strcat(name,inp->name);
	if(!(inp->fil=fopen(name,"rb")))
	{
#ifndef _WIN32
	    if(inp->pipe)
	    {
		int count=5; // wait 5 secs for pipe to open
		while(errno==ENOENT && count--)
		{
		    sleep(1); 
		    if((inp->files[i]=fopen(name,"rb")))goto open2; 
		}
	    }
#endif  
	    print_error(FERROR,sname,"Can't open input file",name);
	open2:;
	}
	if(inp->offset>0)
	{inp->pos=-1;skip_fix(inp,0);} 
	}
    }
    inp->end_pos=irint((inp->f_end-inp->f_start)*inp->freq);
    inp->pos=inp->edf_len=0;// start of file
    return 0;
}



void swap2(short *b, int m, int k)
{ swab(b,b,m*k);}
void swap4(short *b, int m, int k)
{
    // reverse byte order for 4 byte integers
    int i,n=m*k;
    swab(b,b,n);
    n/=2;
    for(i=0;i<n;i+=2)
    { 
	register short x;
	x=b[i];
	b[i]=b[i+1];
	b[i+1]=x;
    }
}
  
void swap8(short *b, int m, int k)
{
    // reverse byte order for 4 byte integers
    int i,n=m*k;
    swab(b,b,n);
    n/=2;
    for(i=0;i<n;i+=4)
    { 
	register short x;
	x=b[i];
	b[i]=b[i+3];
	b[i+3]=x;
	x=b[i+1];
	b[i+1]=b[i+2];
	b[i+2]=x;
    }
}

void write_alarm(short alarm,int i)
{
    // it works using only first inp and out
    char time[30];
    double tod;
    time_t t;
    struct tm *tim;
    tod=in[0].f_start+i/in[0].freq;
    t=(time_t)(tod);
    tim=gmtime(&t);

    strftime(time,30,"%a %e %b %Y %H:%M",tim); 
    tim->tm_hour=tim->tm_min=tim->tm_sec=0;
    tod-=mktime(tim)-timezone;

    if(fprintf(ou[0].alarm,"%s\t0.0\t%.4lf\t%hd\n",time,tod,alarm)<1)
    {
	print_error(ERROR,"write_alarm","Writing to alarm file",NULL);
 
    }
}  




  



int sync_cur;

#define HB  0xf0
#define LB  0x0f

void decode_tf(short *buffer, int n, int m)
{
    // decodes tf packed int12 data in place
    // n is the size of single output buf in bytes(=nchans*2)
    // k is the size of packed data in bytes 
    // m is no of samples to decode
  
    int i,k=3*n/4+4;  //size of packed buffer (alarm +sync are 4 bytes)
    BYTE *buf=(BYTE*)buffer+(m-1)*k;
    register ushort *w=(ushort*)buffer + m*n/2+1; // points to last wrd +2 (4 bytes extra)
    short alarm,sync;
    int pos=sync_cur+m-1;
  
    for(i=m;i>0;i--,buf-=k,pos--)
    {
	register BYTE x,*b=buf+k-3; //points to last byte triple
	sync=buf[0];
	if(pos%200!=sync)
	    print_error(WARNING,"decode_tf","sync error ",in[0].name);
	alarm=buf[2];
	if(alarm && ou[0].alarm)write_alarm(alarm,pos);

	while(b>buf+4)
	{
	    *w-- = ((ushort)*(b+1)<<4) | (ushort)(*(b+2)&LB);
	    *w-- = ((ushort)*b << 4) | ((ushort)(*(b+2)&HB)>>4);
	    b-=3;
	}
	//last separately to avoid overwrite
	// b == buf+4 here
	x=*(b+2);
	*w-- = ((ushort)*(b+1)<<4) | (ushort)(x&LB); // writes in b+2 and b+3
	*w-- = ((ushort)*b << 4) | ((ushort)(x&HB)>>4);// writes in b and b+1
    }
    // shift buffer by 4 bytes (alarm +sync)
  
    memmove(buffer,buffer+2,n*m);
    // last sync and alarm saved at end but overwritten by next read
    buffer[m*n/2]=sync;  // save sync for last sample
    buffer[m*n/2+1]=alarm;
} 
   
void decode_tf_swap(short *buf,int k, int m)
{ int n=k*m;
 decode_tf(buf,k,m);
 swab(buf,buf,n+4);
}

void compute_avg(struct inpdef *inp,double *db,int nn)
{
    struct mont *mtg=&inp->mtg;
    int ii,i,j,n,k,mns,m=inp->no_of_vars,ns=mtg->mchans+1;
    mns=m*(ns+mtg->navg);
    for(ii=0;ii<nn;ii++)
    {
	//for avg0  no lst needed
	if((n=mtg->avg_no[0]))
	    for(k=0;k<m;k++)
	    {
		double sum=0.l;
		for(j=1;j<=n;j++)
		    sum+=db[j*ns+k];
		db[ns*m+k]=sum/n;
	    }
	for(i=1;i<mtg->navg;i++)
	{ 
	    int *lst=mtg->lst_avg[i];
	    if((n=mtg->avg_no[i]))
		for(k=0;k<m;k++)
		{
		    double sum=0.l;
		    for(j=0;j<n;j++)
			sum+=db[(*lst++)*m+k];
		    db[(ns+i)*m+k]=sum/n;
		}
	}
	db+=mns;
    }	
}


void select_chans_short_calib(struct inpdef *inp,char*sbuf,double* db,int n)
{
    // puts selected chans in double buffer dbuf for n points
    // calibrates, assumes no_of_vars==1
    int i,j,m=0;
    short *buf=(short*)sbuf;
    double *dbuf=db;
    struct mont *mtg=&inp->mtg;
    for(j=0;j<n;j++)
    { 
	if(mtg->type>MTG_SEL) 
	    *dbuf++=0.0l;
	for(i=0;i<mtg->mchans;i++)
	{ register int l=mtg->ch_no[i]+m;
	*dbuf++=buf[l]*mtg->chan_calib[i+1];
	} 
	dbuf+=mtg->navg; 
	m+=inp->nchans;
    }
    if(mtg->navg)compute_avg(inp,db,n); 
}

void select_chans_ushort_calib(struct inpdef *inp,char* sbuf,double* db,int n)
{
    // puts selected chans in double buffer dbuf for npoints
    // calibrates, assumes no_of_vars==1
    int i,j,m=0;
    unsigned short *buf=(unsigned short*)sbuf;
    double *dbuf=db;
    struct mont *mtg=&inp->mtg;
    for(j=0;j<n;j++)
    {
	if(mtg->type>MTG_SEL) 
	    *dbuf++=0.0l;
	for(i=0;i<mtg->mchans;i++)
	{ register int l=mtg->ch_no[i]+m;
	*dbuf++=buf[l]*mtg->chan_calib[i+1];
	} 
	dbuf+=mtg->navg;
	m+=inp->nchans;
    }
    if(mtg->navg)compute_avg(inp,db,n); 
}

void get_vars(struct inpdef *inp, void *buf,int n, int ch, int var, int nvars)
{
  
    // gets nvars starting at var in channel ch in point n from rbuf
    char *s;
    if(inp->mtg.type>MTG_SEL)
    {
	print_error(ERROR,"get_vars","can't be used with ref montages",NULL);
	 
    }
    s=(char*)inp->rbuf+inp->bytes_per_var*var
	+inp->bytes_per_chan*(inp->mtg.ch_no[inp->mtg.chan_no[ch]-1]+n*inp->nchans);
    memcpy(buf,s,nvars*inp->bytes_per_var);
}  

// macro SELECT_CHANS expects proper definition of input buf
// no calibration
// puts selected chans in double buffer dbuf
#define SELECT_CHANS { \
  if(inp->dbuf==(double*)inp->rbuf)return;{\
  register int j,m,jj; \
  double *dbuf=db; \
  struct mont *mtg;\
  mtg=&inp->mtg;   \
  jj=mtg->navg*inp->no_of_vars;\
  for(j=0,m=0;j<n;j++) \
  { register int i;\
if(mtg->type>MTG_SEL)  \
  { memset(dbuf,0,inp->no_of_vars*sizeof(double)); \
dbuf+=inp->no_of_vars;}\
for(i=0;i<mtg->mchans;i++) \
{ register int k,l;\
  l=(mtg->ch_no[i]+m)*inp->no_of_vars; \
  for(k=0;k<inp->no_of_vars;k++)   \
	*dbuf++=buf[l+k];  \
}  \
dbuf+=jj;  \
m+=inp->nchans;\
  }\
  if(mtg->navg)compute_avg(inp,db,n);  }   \
}


/*
  void select_chans(struct inpdef *inp, char* sbuf,double *db,int n)
  {
  // special case NO_MTG and fmt is double  no conversion needed
  // no zeros for db 
  // if db= NULL on input just copy of pointer
  if(!db)
  db=(double*)sbuf;
  else
  memcpy(db,sbuf,inp->bytes_win_size*n);
  }
*/
void select_chans_byte(struct inpdef *inp,char* sbuf,double*db,int n)
{
    unsigned char *buf=(unsigned char*)sbuf;
    SELECT_CHANS;
}

void select_chans_short(struct inpdef *inp,char*sbuf,double *db,int n)
{
    short *buf=(short*)sbuf;
    SELECT_CHANS;
}
void select_chans_int(struct inpdef *inp,char*sbuf,double *db,int n)
{
    int *buf=(int*)sbuf;
    SELECT_CHANS;
} 

void select_chans_char(struct inpdef *inp,char* sbuf,double *db,int n)
{
    char *buf=sbuf;
    SELECT_CHANS;
}

void select_chans_ushort(struct inpdef *inp,char*sbuf,double *db,int n)
{
    unsigned short *buf=(unsigned short*)sbuf;
    SELECT_CHANS;
}
void select_chans_uint(struct inpdef *inp,char*sbuf,double *db,int n)
{
    unsigned int *buf=(unsigned int*)sbuf;
    SELECT_CHANS;
} 
void select_chans_float(struct inpdef *inp,char*sbuf,double *db,int n)
{
    float *buf=(float*)sbuf;
    SELECT_CHANS;
} 

void select_chans_double(struct inpdef *inp,char*sbuf,double *db,int n)
{
    double *buf=(double*)sbuf;
 
    SELECT_CHANS;
} 

 
int read_next_win(struct inpdef *inp,int npoints,int part)
{
    // reads next window of data into inp->rbuf
    // selects channels into dbuf
    // works for fixed length (for var_len & npoints=1 it also works)
    // returns no of windows(points) read
    // shifts buffer if ov overlap
    // updates f_cur = begining of current window
    // pos= current position in file in no of windows
    // if part=1 reads partial data, otherwise ignore partial and reads next
    // if eof on stream inp->eof=1; and 0 read
    // sets new_out=1 if first read in new epoch or after gap
 
    int np,n,k;
    char *sbuf;
    double *dbuf;
    int ovpoints;
    int pos;
    int flag=(inp->dbuf && inp->dbuf != (double*)inp->rbuf);
    char *sname="read_next_window";

    if(npoints<=0 )
    {
	char err[50];
  
	npoints=inp->npoints;
	// sprintf(err,"npoints changed to %d npoints",npoints);
	// print_error(WARNING,sname,err,NULL);
    }
    else if(npoints>inp->npoints)
    { char err[50];
  
    new_inp_size(inp,npoints);
    sprintf(err,"Input buffer size increased to %d points",npoints);
    print_error(WARNING,sname,err,NULL);
    }
  
 
    inp->eof=0;
  
    if(inp->gap)
    {
	// start new bufs and out file
	inp->cur=inp->pos;
    }
    else
    { 
	inp->cur+=inp->spoints; 
    }

  
    inp->new_out=inp->gap;

    ovpoints=inp->pos-inp->cur;

    if(ovpoints>0)
    {
	sbuf=(char*)inp->rbuf+inp->spoints*inp->bytes_win_size;
	memmove(inp->rbuf,sbuf,ovpoints*inp->bytes_win_size);
	sbuf =(char*)inp->rbuf+ovpoints*inp->bytes_win_size;
	if(flag)
	{
	    dbuf=inp->dbuf+inp->spoints*inp->dwin_size;
	    memmove(inp->dbuf,dbuf,ovpoints*inp->dwin_size*sizeof(double));
	    dbuf=inp->dbuf+ovpoints*inp->dwin_size;
	}
	else
	    dbuf=NULL;
	np=inp->spoints;
	k=npoints-np;
	pos=inp->pos; 
    }
    else
    { 
	sbuf=(char*)inp->rbuf;
	dbuf=inp->dbuf;
	np=npoints;
	k=0;
	pos=inp->cur;
    }

    inp->f_cur=inp->f_start+inp->cur/inp->freq;


    do
    { 
  
	if(find_next(inp,pos)) // check if not out of epoch sets inp->gap, inp->new_epoch
	{
	    inp->eof=1;
	    if(!part)
		return 0;
	    else
		return k;
	}
	if(inp->gap)
	{
	    // not continue or new epoch
	    if(part)
	    {
		return k;  //it may return 0
	    }
	    else
	    {
		// ignore partial read start new
		sbuf=(char*)inp->rbuf;
		dbuf=inp->dbuf;
		np=npoints;
		k=0;
		inp->cur=inp->pos;
	 
		inp->f_cur=inp->f_start+inp->cur/inp->freq;
	  
		inp->new_out=1;
		inp->gap=0;  // for not partial read gap is always 0 on out
                inp->sel_out_id=inp->sel_id;
                inp->sel_out_text=inp->sel_text;
	    }
	}

	sync_cur=inp->pos; // this is needed for int12
	if((n=inp->read_fun(inp,sbuf,np))>0) // it should be always >0 
	{
	    if(inp->pos+n>=inp->e_end_pos)
		n=inp->e_end_pos-inp->pos;
	    inp->pos+=n;
	    k+=n;
	    if(inp->swap_fun!=NULL)
		inp->swap_fun((short*)sbuf,inp->bytes_win_size,n);
	    if(flag)
	    {
		inp->select_fun(inp,sbuf,dbuf,n);
		dbuf+=n*inp->dwin_size;
	    }
	    sbuf+=n*inp->bytes_win_size;	 	 
	}
	else
	{
	    print_error(WARNING,"read_next_win","unexpected eof",inp->name);
	    inp->eof=1;
	    if(!part)
		return 0;
	    else
		return k;
	}
	  
	np-=n; 
	pos=inp->pos;  
    }while(np>0);
    inp->end=inp->cur+k;  
    return k;
}



void write_chan_data(struct outdef *out,void *buf,int n,int chan,int nvars )
{
    // writes data at current position of out file
    // n number of vars to write from buf in this call
    // if nvars ==0 only write no checks
    // if nvars nvars=<0 assumes nvars=no_of_vars
    // for whole multichannel buffer size= no_of_vars*no_of_chans*npoints*size of type of data
    // if chan==0 && nvars it checks for max size of file and updates count 
    // (convenient for channel by channel analysis starting from 0
 
    // out->new_out is ==1 on return if opened new file (-1 if due to max_len)
    // if in[out->inp_no]->new_out==1 closes last and opens new
    FILE* fil; 
 
    if(n<1)return;
  
    if(nvars)
    {
	if(chan==0)
	    out->new_out=0;
	if(!parallel_read)
	    out->inp_no=cur_inp;
	if (((out->inp_no>-1)&& in[out->inp_no].new_out) || out->bytes >= out->max_len ||
	    (!out->files && !out->fil) || (out->files && !out->files[chan]))
	{
	    if(out->bytes>=out->max_len && out->offset>0)
	        print_error(FERROR,"write_chan_data","written bytes>Max_len and offset>0.Increase Max_len for output.",out->name);
	    if(out->files)
		new_out_file_chan(out,chan);
	    else if(chan==0)
		new_out_file(out);
	}

	if (out->bytes >= out->max_len)
	    out->new_out=-out->new_out;
    }
 
    if(nvars<0)
	nvars=out->no_of_vars;
    if(out->files)
	fil=out->files[chan];
    else
	fil=out->fil;

    if(nvars && out->var_len && (chan==0 ||out->files))
    {
	if(out->swap_fun)
	    swap4((short*)&nvars,4,1);
	fwrite(&nvars,4,1,fil);
    }

    if(out->swap_fun)
	out->swap_fun(buf,out->bytes_per_var,n);
    if((int)fwrite(buf,out->bytes_per_var,n,fil)<n)
    {
	print_error(FERROR,"write_chan_data","writing output file",out->name);
  
    }
    if(out->pipe)fflush(fil);  
 
    if(nvars && chan==0)
	out->count++;

    if(!out->files)
	out->bytes+=n*out->bytes_per_var;
    else
    {
	int i;
	out->mbytes[chan]+=n*out->bytes_per_var;
	if(chan==(out->no_of_chans-1))
	    for(i=0;i<out->no_of_chans;i++)
		if(out->bytes<out->mbytes[i])
		    out->bytes=out->mbytes[i];
    }

}

void write_data_raw(struct outdef *out, void *buf, int npoints, int nvars,int varsize)
{ // write no_of_chans starting from start_chan
  // nvars no of vars per channel if nvars <1 nvars=out->no_of_vars
  // it assumes all chans have same nvars
    // varsize is size of data in buffer i.e. varsize*out->no_o_vars skips to next chan
    char *b,*bb;
    int i,j,nv,nc,n;
    //  struct inpdef *inp=in+out->inp_no;
    if(npoints<1)return;
    if(nvars<1)
	nvars=n=out->no_of_vars;  //this is default  
    else
	out->no_of_vars=n=nvars; 
    if(out->var_len)
	npoints=1;
    if(varsize<1)varsize=out->bytes_per_var;
    nv=varsize*nvars; 
    nc=nvars*out->all_chans*varsize; 
    bb=(char*)buf+out->start_chan*nvars*varsize;
    for(j=0;j<npoints;j++,n=0,bb+=nc)	
    {   
	for(i=0,b=bb;i<out->no_of_chans;i++,b+=nv)
	{ if(out->no_chan_vars)   
	    write_chan_data(out,b,out->no_chan_vars[i],i,out->no_chan_vars[i]); 
	else
	    write_chan_data(out,b,nvars,i,n); 
	}
	if(!n)out->count++; 
    }  
}

void close_io()
{
    // writes last list file line(if not pipe)
    // closes all open files
  
    int i,j;
    char name[MAX_FILE_NAME],*sname="close_io";
 
    for(i=0;i<no_of_inp;i++)
    {
	if(in[i].files)
	    for(j=0;j<in[i].mtg.mchans;j++)
	    {
		if(in[i].files[j] && fclose(in[i].files[j]))
		    print_error(WARNING,sname,"error closing inp file",in[i].name);
	 
	    }
	else
	    if(in[i].fil && fclose(in[i].fil))
		print_error(WARNING,sname,"error closing inp file",in[i].name);
	 

	if(in[i].pipe)
	{
	    strcpy(name,in[i].path);
	    strcat(name,in[i].name);
	    if(unlink(name))  // delete pipe
		print_error(WARNING,sname,"Can't unlink inp pipe",name);
	}
	if(in[i].lst && fclose(in[i].lst))
	    print_error(WARNING,sname,"error closing inp list file",in[i].lst_name);
	
	if(in[i].pipe)
	{ 
	    strcpy(name,in[i].path);
	    strcat(name,in[i].lst_name);
	    if(unlink(name))  // delete pipe
		print_error(WARNING,sname,"Can't unlink list pipe",name);
	}
	if(in[i].sel)
	{
	    if(fclose(in[i].sel))
		print_error(WARNING,sname,"error closing select file",in[i].sel_name);
	}
    }
  
    for(i=0;i<no_of_out;i++)
    {
	if(ou[i].lst && !(ou[i].pipe))
	{
	    if(ou[i].rate_flag)
		ou[i].len=(ou[i].count-1)/ou[i].win_rate+ou[i].win_len;
	    else
		ou[i].len=(ou[i].count-1)*ou[i].win_shift+ou[i].win_len;
	 if(ou[i].inp_no>-1)
          {
            ou[i].sel_id=in[ou[i].inp_no].sel_id;
	    ou[i].sel_text=in[ou[i].inp_no].sel_text;
          }
	    write_list_line(&ou[i]);
	}
  
	if(ou[i].files)
	    for(j=0;j<ou[i].no_of_chans;j++)
	    {
		if(ou[i].files[j] && fclose(ou[i].files[j]))
		    print_error(WARNING,sname,"error closing out file",ou[i].name);
	    }
	else
	    if(ou[i].fil && fclose(ou[i].fil))
		print_error(WARNING,sname,"error closing out file",ou[i].name);

	if(ou[i].lst && fclose(ou[i].lst))	
	    print_error(WARNING,sname,"error closing out list file",ou[i].lst_name);
  
    }
    if(ou[0].alarm)fclose(ou[0].alarm);

    print_error(INFO, "main", "Done", NULL);
}




int find_next(struct inpdef *inp,int n)
{ 
    // replacement for find_forward in  read_next_win
    // works on epochs(from list file) in memory
    // alows for overlapping selections and files
    // return -1 when eof 0 otherwise
    // positions input at n or begining of next epoch
    // sets gap, new_epoch
    // calls pos_file
 

    int i,k,m=2*inp->list[0],eflag=1,nlen;
    double *list=inp->list+1,sec,sec_end,last,eps=0.99/inp->samp_rate;
//    printf("n=%d\n",n);
    if(inp->pipe)
    {
	if(n>=inp->pos)
	    return find_forward(inp,n);
	else
	    print_error(FERROR,"find_next","no back for pipe",inp->name);
    }


    inp->gap=inp->new_epoch=0;
    nlen=inp->inp_npoints/inp->inp_spoints-1;
    last=sec=inp->f_start+n/inp->freq;
    sec_end=sec+nlen/inp->freq;
    k=2*inp->cur_epoch-2;

    // check if in epoch

    if(n+nlen>=inp->e_end_pos)
    {
    next_epoch:
	if(read_select_line(inp))return -1;
	if( inp->pos!=inp->e_start_pos || fabs(sec-inp->e_start)>eps )
	{
	    sec=inp->e_start; // set to the begining of next epoch
	    sec_end=sec+nlen/inp->freq;
	    inp->gap=1; // always gap even when continous
	    k=0; //start looking for file from the begining
	    eflag=0;
//            printf("sec=%lf,sec_end=%lf,n=%d,nlen=%d\n",sec,sec_end,n,nlen);
	}  
    }  
    else if(n<inp->e_start_pos)  // just to be sure
{ //           printf("start-pos=%d,sec=%lf,sec_end=%lf,n=%d,nlen=%d\n",inp->e_start_pos,sec,sec_end,n,nlen);
	print_error(FERROR,"find_next","selections out of order",inp->name);
}

    // here in epoch
    // check if in file and go to next if not

    for(i=k;i<m;i+=2)
    {
	
	if((sec-list[i])<=-eps)  // sec<list[i]
	{
	    sec=list[i];
	    sec_end=sec+nlen/inp->freq;
	    if((sec_end-inp->e_end)>-eps)goto next_epoch; // sec_end >= inp->e_end
	    else goto in_file;
	}
	else
	{
	    if((sec_end-list[i+1])<=-eps)goto in_file; // sec_end < list[i+1]
	    else if(eflag && i<m-2)
	    {
		if(sec < list[i+2] || fabs(sec-inp->e_start)>eps )// set to beg of file only if not beg of epoch
		    sec=list[i+2];
		sec_end=sec+nlen/inp->freq;
		if((sec_end-inp->e_end)>-eps)goto next_epoch; // sec_end >= inp_e_end
	    }
	}
    }
    return -1;  // eof

 in_file:
    k=i/2+1; // curent file in list
    inp->gap=inp->gap ||(fabs(last-sec)>eps);
    if(!inp->sel_comment)
    {
	inp->sel_id=inp->id[k];
	inp->sel_text=inp->text[k];
    }
    if(pos_file(inp,k,sec))
	return 0; // OK
    else
	print_error(FERROR,"find_next","pos_file returns 0",inp->sel_name);
 
}
/*
  int find(struct inpdef *inp,double sec)
  { 
  return find_next(inp,irint((sec-inp->f_start)*inp->freq));
  }
*/ 
  


int find_forward(struct inpdef *inp, int n)
{ // finds and opens file containing data for pos n=no_of windows from start
  // returns 0 if OK
  // if eof on list or select returns -1
  // checks if in epoch
  // repositions list file and input file
  // if pos between epochs or before first epoch
  // starts at the next closest epoch
  // it assumes that list file(and select) is open before
  // before calling inp->pos and inp->f_start should be set 
  // positions are measured in windows from f_start of current file
  // assumes n >=0 
  // before first call inp file should be open
  // sets inp->gap and inp->pos (in skip)
  

    inp->gap=inp->new_epoch=0;
    for(;;)
    {
	while(n>=inp->e_end_pos)
	    if(read_select_line(inp))return -1;  // sets new_epoch
	if(n<inp->e_start_pos)
	{
	    n=inp->e_start_pos; // set to the begining of next epoch
	    inp->gap=1;
	}
	if(n<0)
	{
	    // here only after new_file
	    inp->gap=1;
	    if(inp->e_end_pos<=0)
	    {
		if(read_select_line(inp))
		    return -1;  // sets new_epoch
		else
		{
		    n=0;
		    continue;
		}
	    }
	    else	
	    {  
		return 0;
	    }
	}
	// here n in epoch && >0
	if(inp->pipe)
	    if(!inp->skip_fun(inp,n))
		return 0;
	    else
		do
		{
		    double last=inp->f_start;
		    if(new_inp_file(inp,1))return -1;	 
	   
		    n+=irint((last-inp->f_start)*inp->freq);
	   
		    // relative to new f_start
		}
		while(inp->skip_fun(inp,n)); 
	else
	{
	    // no pipe file	  
	    if(n<inp->end_pos) // in epoch and in  file	
		return inp->skip_fun(inp,n);  // return at pos n>=0 (gap=0 or 1)
	    // next file
	    do
	    {
		double last=inp->f_start;
		if(read_list_line(inp))return -1;
	  
		n+=irint((last-inp->f_start)*inp->freq);
	  
		// relative to new f_start
	 
	    } while(n>=inp->end_pos);
	    new_inp_file(inp,0); // opens current; no read from list	   
	}
	 
	// in file but  maybe out of epoch
    }
}
short *edf_tran=NULL;
void transpose(short *x,int m,int n)
{
    // transpose short matrix x(m,n)(row-wise. i.e. x11, x12..)
  
    short *b;
    int i,j;
    memcpy(edf_tran,x,m*n*sizeof(short));
    for(i=0;i<m;i++)
    {
	b=edf_tran+i;
	for(j=0;j<n;j++,b+=m)
	    *x++=*b;
    }
}
  
  
  
int read_win_edf(struct inpdef *inp, char *sbuf, int np)
{
    // provides n input samples from edf file
    int n=np,k=0;
    int rec_len=inp->rec_len; 
    while(n>0)
    {
 
	if(inp->edf_len>0)
	{
	    if(inp->edf_len>n)
	    { int l=inp->bytes_win_size*n;
	    memcpy(sbuf,inp->edf_sbuf,l);
	    inp->edf_len-=n;
	    inp->edf_sbuf+=l/2;
	    n=0;
	    break;
	    }
	    else
	    { int l=inp->bytes_win_size*inp->edf_len;
	    memcpy(sbuf,inp->edf_sbuf,l);
	    n-=inp->edf_len;
	    sbuf+=l;
	    }
	}
	if(inp->zeros)
		k=put_zeros(inp,inp->edf_buf,inp->bytes_win_size,rec_len);
	else
		k=fread(inp->edf_buf,inp->bytes_win_size,rec_len,inp->fil);

	if(k!=rec_len && ferror(inp->fil))
	    print_error(FERROR,"read_win_edf","read error",inp->name);
	else
	{
	    inp->edf_len=k;
	    inp->edf_sbuf=inp->edf_buf;
	    if(k==0) //eof
		break;
	    else	  
		transpose(inp->edf_buf,k,inp->nchans);
	}
    }
    return np-n;
}
	  
	  
  
  


int read_win_fix(struct  inpdef *inp, char* sbuf, int n)
{
    // reads n input windows of fixed length from  one file
    // returns no of windows read
    int k;
    if(inp->zeros)
		k=put_zeros(inp,sbuf,inp->bytes_per_win,n);
	else
		k=fread(sbuf,inp->bytes_per_win,n,inp->fil);
    if(k!=n && ferror(inp->fil))
    {
	print_error(FERROR,"read_win_fix","read error",inp->name);
    }
    return k; 
}

int read_win_fix_chan(struct  inpdef *inp, char*sbuf, int n,int chan)
{
    // reads n input windows of fixed length from  one channel file
    // returns no of windows read
    int k;
    char *buf=sbuf+inp->bytes_per_chan*chan;
    FILE* fil=inp->files[chan];
     if(inp->zeros)
		k=put_zeros(inp,buf,inp->bytes_per_chan,n);
	else
		k=fread(buf,inp->bytes_per_chan,n,inp->fil);
    if(k!=n && ferror(fil))
    {
	sprintf(buf_line,"%s.%3.3d",inp->name,chan);
	print_error(FERROR,"read_win_fix_chan","read error",buf_line);
    }
    return k; 
}
 

int read_win_var(struct  inpdef *inp, char* sbuf, int n)
{
    // reads one window  of variable length from  one file 
    // n is ignored
    // returns 1, always one window read
    // updates no of vars read
    // if too large for buffer trunctates
    char *sname="read_win_var";
    n=0;
    if(inp->zeros)
    {
	  inp->no_of_vars=inp->var_len;
	  put_zeros(inp,sbuf,inp->bytes_per_var,inp->var_len);
	  return 1;
    }
    if(fread(&n,4,1,inp->fil)!=1 && feof(inp->fil))return 0; //eof 
 
    if(inp->swap_fun)swap4((short*)&n,4,1);
    if(n>inp->var_len)
    {  n=inp->var_len;
    print_error(WARNING,sname,"record length trunctated",inp->name);
  
    }
  
    if(n<1 ||(((int)fread(sbuf,inp->bytes_per_var,n,inp->fil))!=n))
    { // eof here is error
	print_error(FERROR,sname,"read error",inp->name);
    }
    else
	inp->no_of_vars=n;
    return 1;
}  
 
int read_win_var_chan(struct  inpdef *inp, char* sbuf, int n,int chan)
{
    // reads one input window of variable length from  one channel
    // updates no of vars read
    // only one is read regardless of n
    FILE* fil=inp->files[chan];
    char *buf=sbuf+inp->bytes_per_chan*chan,*sname="read_win_var_chan";
    n=0;
     if(inp->zeros)
    {
	  inp->no_chan_vars[chan]=inp->var_len;
	  put_zeros(inp,sbuf,inp->bytes_per_var,inp->var_len);
	  return 1;
    }
    if(fread(&n,4,1,fil)!=1 && feof(fil))return 0;  //eof 
    if(inp->swap_fun)swap4((short*)&n,4,1);
    if(n>inp->var_len)
    {  n=inp->var_len;
    print_error(WARNING,sname,"record length trunctated",inp->name);
    }
    if(n<1 ||(((int)fread(buf,inp->bytes_per_var,n,fil))!=n))
    {
	sprintf(buf_line,"%s.%3.3d",inp->name,chan);
	print_error(FERROR,sname,"read error",buf_line);
    }
    else
	inp->no_chan_vars[chan]=n;
    return 1;   
}

int read_win_fix_many(struct inpdef *inp,char* sbuf,int n)
{
    // reads n input windows of fixed length from  all one channel files
    // returns no of windows read  
    // assumes eof in all files symultanously
    int i,k;
    char *s=sbuf;
    if(inp->zeros)
   {
	k=put_zeros(inp,sbuf,inp->bytes_per_win,n);
	return k;
   }
	
    for(i=0;i<n;i++)
	for(k=0;k<inp->nchans;k++,s+=inp->bytes_per_chan)
	    if(!fread(s,inp->bytes_per_chan,1,inp->files[k]))
		goto end;
 end:
    if(i<n && ferror(inp->files[k]))
    {
	sprintf(buf_line,"%s.%3.3d",inp->name,k);
	print_error(FERROR,"read_win_fix_many","read error",buf_line);
    }
    return i;
}

int read_win_var_many(struct inpdef *inp,char* sbuf,int n)
{
    // reads n input windows of var length from  all one channel files
    // returns no of windows read (for var only 1 is supported)
    // assumes eof in all files symultanously
    int i,k;
    for(k=0;k<inp->nchans;k++)
	i=read_win_var_chan(inp,sbuf,1,k);
    return i;
}





int skip_fix(struct inpdef *inp,int n)
{
    // always forward (n>0)
    // skipps current file to position n=no of windows from f_start
    // returns -1 if eof 0 if OK
    // sets inp->pos
    int i;
    long offset;
//    if(n==inp->pos)return 0;
    if(inp->zeros){inp->pos=n;return 0;}
    if(inp->pipe)
    {
	if(n>0)
	{
	    for(i=0;i<n;i++)
		if(read_win_fix(inp,inp->rbuf,1)!=1)return -1; // eof
	    inp->pos=n;
	}
	return 0;
    }

    // file
    offset=n*inp->bytes_per_win+inp->offset;
    if(fseek(inp->fil,offset,SEEK_SET))
    { 
	print_error(FERROR,"skip_fix","fseek error",inp->name);
    }
    inp->pos=n;
    return 0; // OK
}
  
int skip_edf(struct inpdef *inp,int n)
{
    // always forward (n>0)
    // skipps current edf file to position n=no of windows from f_start
    // returns -1 if eof 0 if OK
    // sets inp->pos
    int i,rec_len,nrec;
    long offset;
//    if(n==inp->pos)return 0;
    if(inp->zeros){inp->pos=n;return 0;}
    if(inp->pipe)
    {
	if(n>0)
	{
	    for(i=0;i<n;i++)
		if(read_win_edf(inp,(char*)edf_tran,1)!=1)return -1; // eof
	    inp->pos=n;
	    inp->edf_len=0;
	}
	return 0;
    }
    // file
    rec_len=inp->rec_len;
    nrec=n/rec_len;
    offset=nrec*rec_len*inp->bytes_win_size+inp->offset;

    if(fseek(inp->fil,offset,SEEK_SET))
    { 
	print_error(FERROR,"skip_edf","fseek error",inp->name);
    }
    inp->edf_len=0;
    if(i=n-nrec*rec_len)
	if(read_win_edf(inp,(char*)edf_tran,i)!=i)return -1; // eof
    inp->pos=n;
    return 0; // OK
}


 
int skip_fix_many(struct inpdef *inp,int n)
{
    int i,k;
    long offset;

//    if(n==inp->pos)return 0;
    if(inp->zeros){inp->pos=n;return 0;}
    if(inp->pipe)
    {   
	if(n>0)
	{
	    for(i=0;i<n;i++)
		for(k=0;k<inp->nchans;k++)
		    if(read_win_fix_chan(inp,inp->rbuf,1,k)!=1)return -1; // eof
	    inp->pos=n;
	}
	return 0;
    } 

    // files
    offset=n*inp->bytes_per_chan+inp->offset;
    for(k=0;k<inp->nchans;k++)
	if(fseek(inp->files[k],offset,SEEK_SET))
	{
	    sprintf(buf_line,"%s.%3.3d",inp->name,inp->fil_ch_no+inp->ch_no[k]);
	    print_error(FERROR,"skip_fix_many","fseek error",buf_line);
	
	}
    inp->pos=n;
    return 0; // OK
}

int skip_var_many(struct inpdef *inp,int n)
{ // is no of windows from start of file  
    int i,k;
//    if(n==inp->pos)return 0;
    if(inp->zeros){inp->pos=n;return 0;}
    n-=inp->pos;
    for(i=0;i<n;i++)
	for(k=0;k<inp->nchans;k++)
	    if(read_win_var_chan(inp,inp->rbuf,1,k)!=1)return -1; // eof
    inp->pos+=n;
    return 0;
}

int skip_var(struct inpdef *inp,int n)
{
    int i;
//    if(n==inp->pos)return 0;
    if(inp->zeros){inp->pos=n;return 0;}
    n-=inp->pos;
    for(i=0;i<n;i++)
	if(read_win_var(inp,inp->rbuf,1)!=1)return -1; // eof
    inp->pos+=n;
    return 0;
}






void create_all_labels(int last_inp)
{
    int i,j,k,flag=0;
    char *sname="create_all_labels";
			 
    all_labels=calloc(all_chans+1,sizeof(char*));
    if(!all_labels)mem_error(sname,"all_labels");
    set_lab_no(all_labels,all_chans);
    all_units=calloc(all_chans+1,sizeof(char*));
    if(!all_units)mem_error(sname,"all_units");
    set_lab_no(all_units,all_chans);
  
    if(last_inp<1 )
    {
	//create default labels
	for(j=1;j<all_chans;j++)
	{   
	    char buf[100];
	    sprintf(buf,"Chan_%d",j);
	    all_labels[j]=strdup(buf);
	}
    }
    else
    {

	// concatenate labels and units from all previous inputs (needed for join) 
			 
	for(i=0,k=1;i<last_inp;i++)
	    for(j=1;j<=in[i].no_of_chans;j++,k++)
	    {
                if(in[i].mtg.chan_labels)
		all_labels[k]=in[i].mtg.chan_labels[j];
	        else
                {char buf[100];
	          sprintf(buf,"Chan_%d",k);
	           all_labels[k]=strdup(buf);
                }
		if(in[i].mtg.chan_units)
		{all_units[k]=in[i].mtg.chan_units[j];flag=1;}
					  
	    }
    }
    if(!flag)
    {free(all_units);all_units=NULL;}
}

int read_cfg(char *fname, struct section *head)
{
    // reads definitions of input output and rest if head not NULL
    // fills global in  & ou struct  or arrays of struct
    // it is called first so no header read yet

    FILE * fil;
    struct section *s;
    char *sname="read_cfg";
    int ret,ni=0,no=0,nn=0,n,last_inp=0;
    char name[MAX_FILE_NAME],buf[MAX_LINE];
 
    fil=fopen(fname,"r");
    if(fil==NULL){
	print_error(FERROR,sname,"Can't open ctl file",fname);
    }
    ctl_name=strdup(fname);
    // general defaults
    no_of_inp=no_of_out=1;
    in=&inp_def;
    ou=&out_def;
    inp_dflt.calib=1;
    inp_ops_fil=fil;
    ret=fscanf(inp_ops_fil,"%99[^\n]\n",buf_line);
    copy_ops=0;
    all_chans=0;
    all_labels=NULL;
    all_units=NULL;

    while(!feof(fil)
	  && (ni <no_of_inp || no<no_of_out))
    {
	if(ret<1 || (ret=skip_sec(name))<1){
	    print_error(FERROR,sname,"read ctl file",fname);

	}
 
	if((s=lookup_sec(&inp_out_sec,name)))
	{
	    if(s==&inp_out_sec)
	    {
		if(nn)
		{
		    char erbuf[100];
		    sprintf(erbuf,"section %s in ctl %s ignored",name,fname);
		    print_error(WARNING,sname,erbuf,NULL);
		    ret=fscanf(inp_ops_fil,"%99[^\n]\n",buf_line);
		    if(ret>0 && copy_ops)fprintf(out_ops_fil,"%s\n",buf_line);
		}
	
	
		else
		{
		    ret=read_sec(s);
		    nn=1;
		    if(no_of_inp>1)
			if(!(in=calloc(no_of_inp,sizeof(struct inpdef))))
			    mem_error(sname,"ctl in");
		    if(no_of_out>1)
			if(!(ou=calloc(no_of_out,sizeof(struct outdef))))
			    mem_error(sname,"ctl ou");
		}
	    }
	    else if(s==&data_inp_sec)
	    {
		nn=1;
		ni++;
		if(no_of_inp>0 && ni>no_of_inp)
		{
		    char erbuf[100];
		    sprintf(erbuf,"%d section %s in ctl %s ignored",ni,name,fname);
		    print_error(WARNING,sname,erbuf,NULL);
		    ret=fscanf(inp_ops_fil,"%99[^\n]\n",buf_line);
		    if(ret>0 && copy_ops)fprintf(out_ops_fil,"%s\n",buf_line);
		}
		else
		{
		    // defaults
		    memcpy(&inp_def,&inp_dflt,sizeof(struct inpdef)); 
		    ret=read_sec(s);
		    last_inp=ni;
		    if(inp_def.hdr_name==NULL)
		    {
			print_error(FERROR,sname,"INPUT header file name is missing",fname);
		  
		    }
		    if(!inp_def.path)
			inp_def.path=strdup("./");
		    else if(inp_def.path[strlen(inp_def.path)-1]!='/')
		    {
			char *tmp=malloc(strlen(inp_def.path)+2);
			strcpy(tmp,inp_def.path);
			strcat(tmp,"/");
			free(inp_def.path);
			inp_def.path=tmp;
		    }
		
		    // now read header  
		  
		    if(*inp_def.hdr_name=='/'||!strncmp(inp_def.hdr_name,"./",2))
			strcpy(name,inp_def.hdr_name);
		    else
		    {
			strcpy(name,inp_def.path);
			strcat(name,inp_def.hdr_name);
		 
		    }
		    memcpy(buf,buf_line, MAX_LINE);
		    read_hdr(name,hdr_sec_list);
		 
		    memcpy(buf_line,buf,MAX_LINE);
		    inp_ops_fil=fil;
		    inp_def.lst_name=strdup(h_fmt.lst_name);
		    if(!inp_def.title)
		    {
			if(h_fmt.title)inp_def.title=strdup(h_fmt.title);
			else inp_def.title=" ";
		    }
		    if(h_fmt.fil_ch_no>0) inp_def.fil_ch_no=h_fmt.fil_ch_no;
		    if(h_fmt.bytes_per_win<=0)
		    {
			print_error(FERROR,sname,"bytes_per_window <=0 in hdr",name);
		  
		    }
		    inp_def.bytes_per_win=h_fmt.bytes_per_win;

		    if(h_fmt.file_fmt==text)
		    {
			print_error(FERROR,sname,"Text format not implemented",name);
		 
		    }
		    if(h_fmt.file_fmt<0)
		    {
			print_error(FERROR,sname,"Unknown file format in hdr",name);
		 
		    }
		    if(h_fmt.byte_order<0)
		    {
			print_error(FERROR,sname,"Unknown byte order in hdr",name);
		  
		    }

		    if(h_fmt.nchans<1)
		    {
			print_error(FERROR,sname,"no of chans <1 in hdr",name);
		 
		    }

		    inp_def.nchans=h_fmt.nchans;

		    if(h_fmt.win_len<=0.)
		    {
			print_error(FERROR,sname,"win_len <1 in hdr",name);
		 
		    }  
		
		    if(h_fmt.win_rate>0.)
		    {
			h_fmt.win_shift=1.0l/h_fmt.win_rate;
			inp_def.rate_flag=1;
		    }
		    else if(h_fmt.win_shift>0.)
		    {
			h_fmt.win_rate=1.0l/h_fmt.win_shift; 
			inp_def.rate_flag=0;
		    }
			  
		    else
		    {
			print_error(FERROR,sname,"win_shift and win_rate<1 in hdr",name);
		    }

		    if(inp_def.win_len<=0.)
			inp_def.win_len=h_fmt.win_shift; //dflt

		    if(inp_def.no_of_vars==0)
			inp_def.no_of_vars=h_fmt.no_of_vars; //dflt
		  

		    if(inp_def.no_of_vars<0)
		    {
			inp_def.no_of_vars=-inp_def.no_of_vars;
			inp_def.var_len=inp_def.no_of_vars;
		    }
		    else if(inp_def.no_of_vars==0)
		    {
			print_error(FERROR,sname,"no_of_vars==0 in hdr",name);
		 
		    }
		    else
			inp_def.var_len=0;
		  
		
		    if(inp_def.var_len)
		    {
			inp_def.npoints=1; // for variable length only 1 input win is possible
			inp_def.spoints=1;
			print_error(WARNING,sname,"npoints & spoints changed to 1 for var len records",fname);
		 
		    }		  
		    inp_def.inp_npoints=irint(h_fmt.win_len*h_fmt.samp_rate);
		    inp_def.inp_spoints=irint(h_fmt.win_shift*h_fmt.samp_rate);
		 
		 
		  
		    if(inp_def.npoints>0)
		    {
			if(inp_def.rate_flag)
			    inp_def.win_len=inp_def.npoints/h_fmt.win_rate;
			else
			    inp_def.win_len=inp_def.npoints*h_fmt.win_shift;
		    }
		    else if(inp_def.win_len>0 )
		    {
			if(inp_def.rate_flag)
			{
			    inp_def.npoints=irint(inp_def.win_len*h_fmt.win_rate);
			    inp_def.win_len=inp_def.npoints/h_fmt.win_rate;
			}
			else
			{
			    inp_def.npoints=irint(inp_def.win_len/h_fmt.win_shift);
			    inp_def.win_len=inp_def.npoints*h_fmt.win_shift;
			}
		    }

		    else   
		    {
			inp_def.npoints=1;
			inp_def.win_len=h_fmt.win_shift;
		    }
		  
		  
		    if(inp_def.spoints>0)
		    {	
			if(inp_def.rate_flag) 
			{
			    inp_def.win_rate=h_fmt.win_rate/inp_def.spoints;
			    inp_def.win_shift=1.0l/inp_def.win_rate;
			}
			else
			{
			    inp_def.win_shift=inp_def.spoints*h_fmt.win_shift;
			    inp_def.win_rate=1.0l/inp_def.win_shift;
			}	
		    }
		    else if(inp_def.win_shift>0)
		    {
			 
			if(inp_def.rate_flag)
			{
			    inp_def.spoints=irint(inp_def.win_shift*h_fmt.win_rate);
			    inp_def.win_shift=inp_def.spoints/h_fmt.win_rate;
			    inp_def.win_rate=h_fmt.win_rate/inp_def.spoints;
			}
				  
			else
			{
			    inp_def.spoints=irint(inp_def.win_shift/h_fmt.win_shift);
			    inp_def.win_shift=inp_def.spoints*h_fmt.win_shift;
			    inp_def.win_rate=1.0l/inp_def.win_shift;
			}

			 
			 
		    }
		    else 
		    {
			inp_def.spoints=inp_def.npoints;
			inp_def.win_shift=inp_def.win_len;
			inp_def.win_rate=1.0l/inp_def.win_shift;
		    }
		
		    if(inp_def.win_shift==h_fmt.win_shift)
			inp_def.win_rate=h_fmt.win_rate;
		   
//	inp_def.rate_flag=prec(inp_def.win_rate)>prec(inp_def.win_shift);

		 
		    inp_def.freq=h_fmt.win_rate;
		    inp_def.samp_rate=h_fmt.samp_rate;
		    inp_def.offset=h_fmt.offset;
		 
		    inp_def.file_fmt=h_fmt.file_fmt;
		    inp_def.new_out=0;
		    inp_def.f_cur=0.;
		    inp_def.pos=0;
		  
		   
		    if(h_fmt.file_fmt==bin_edf && (h_fmt.nchans==1 || h_fmt.no_of_vars==1)){
			inp_def.file_fmt==bin_short;
			}

		    switch(inp_def.file_fmt)
		    {
			case bin_char:
			case bin_byte:
			    inp_def.bytes_per_var=1;
			    break;
			case bin_int12:
			case bin_short:
			case bin_ushort: 
			case bin_edf:
			    inp_def.bytes_per_var=sizeof(short);
			    break;
			case bin_int:
			case bin_uint:
			    inp_def.bytes_per_var=sizeof(int);
			    break;
			case bin_float:
			    inp_def.bytes_per_var=sizeof(float);
			    break;
			case bin_double:
			    inp_def.bytes_per_var=sizeof(double);
		    }  
		 
		    if(inp_def.mtg_name)
		    {
			memcpy(buf,buf_line, MAX_LINE);
			if(*inp_def.mtg_name=='/'||!strncmp(inp_def.mtg_name,"./",2))
			    strcpy(name,inp_def.mtg_name);
			else
			{
			    strcpy(name,inp_def.path);
			    strcat(name,inp_def.mtg_name);
			}
			read_mtg(name,&inp_def.mtg,&h_fmt);
			memcpy(buf_line,buf,MAX_LINE);
		    }
		    else
			read_mtg(NULL,&inp_def.mtg,&h_fmt);
		    if(inp_def.var_len)
			v_fmt.nchans=1;
		    else
			v_fmt.nchans=inp_def.no_of_vars;

		    if(inp_def.vmtg_name)
		    {
			memcpy(buf,buf_line, MAX_LINE);
			if(*inp_def.vmtg_name=='/'||!strncmp(inp_def.mtg_name,"./",2))
			    strcpy(name,inp_def.vmtg_name);
			else
			{
			    strcpy(name,inp_def.path);
			    strcat(name,inp_def.vmtg_name);
			 			  
			}
			read_mtg(name,&inp_def.vmtg,&v_fmt);
			memcpy(buf_line,buf,MAX_LINE);
		    }
		    else
			read_mtg(NULL,&inp_def.vmtg,&v_fmt);


		    inp_def.no_of_chans=inp_def.mtg.nchans;
		    if(inp_mode==mode_parallel)
			all_chans+=inp_def.no_of_chans; 
                    else
			all_chans=inp_def.no_of_chans; //for serial it is always last
		 
		  
		    if(h_fmt.chans_per_file==-1)
		    {
			int j;
			inp_def.nchans=inp_def.mtg.mchans;  // only needed channels are read
			if(!(inp_def.files=(FILE**)calloc(inp_def.nchans,sizeof(FILE*))))
			    mem_error(sname,"inp files");
			if(!(inp_def.ch_no=(int*)calloc(inp_def.nchans,sizeof(int))))
			    mem_error(sname,"inp files ch_no");
		 
			for(j=0;j<inp_def.nchans;j++)
			{
			    inp_def.ch_no[j]=inp_def.mtg.ch_no[j];
			    inp_def.mtg.ch_no[j]=j;
			}
		    }
		    else
			inp_def.files=NULL;
  
		 
		  
 
		    if(inp_def.var_len)
		    {
			if(inp_def.mtg.type>MTG_SEL || inp_def.vmtg.type>MTG_SEL)
			{
			    print_error(FERROR,sname,"Only selection montage supported for var len data",name);
			}
		   
			if(inp_def.files)
			{
			    inp_def.skip_fun=skip_var_many;
			    inp_def.read_fun=read_win_var_many;
			    if(!(inp_def.no_chan_vars=(int*)calloc(inp_def.no_of_chans,sizeof(int))))
				mem_error(sname,"inp no_chan_vars");
			}
			else
			{
			    inp_def.skip_fun=skip_var;
			    inp_def.read_fun=read_win_var;
			 
			}
		    }
		    else 
		    {
			if(inp_def.files)
			{ 
			    inp_def.skip_fun=skip_fix_many;
			    inp_def.read_fun=read_win_fix_many;
			}
			else if(inp_def.file_fmt==bin_edf)
			{
			    inp_def.skip_fun=skip_edf;
			    inp_def.read_fun=read_win_edf;  	 
			}
			else 
			{
			    inp_def.skip_fun=skip_fix;
			    inp_def.read_fun=read_win_fix;  	 
			}
		   
		    }
		   
		    if(inp_def.file_fmt==bin_int12)
			inp_def.swap_fun=decode_tf;
		    else
			inp_def.swap_fun=NULL;

		    if(h_fmt.byte_order != MY_BYTE_ORDER)
			switch(inp_def.file_fmt)
			{
			    case bin_short:
			    case bin_ushort: 
			    case bin_edf:
				inp_def.swap_fun=swap2;
				break;
			    case bin_int:
			    case bin_uint:
			    case bin_float:
				inp_def.swap_fun=swap4;
				break;
			    case bin_double:
				inp_def.swap_fun=swap8;
			}  
		    switch(inp_def.file_fmt)
		    {
			case bin_char:
			    inp_def.select_fun=select_chans_char;
			    break;
			case bin_byte:
			    inp_def.select_fun=select_chans_byte;
			    break;
			case bin_int12:
			case bin_short:
			case bin_edf:
			    if(inp_def.calib && inp_def.mtg.chan_calib)
				inp_def.select_fun=select_chans_short_calib;
			    else
				inp_def.select_fun=select_chans_short;
			    break;
			case bin_ushort:  
			    if(inp_def.calib && inp_def.mtg.chan_calib)
				inp_def.select_fun=select_chans_ushort_calib;
			    else
				inp_def.select_fun=select_chans_ushort;
			    break;
			case bin_int:
			    inp_def.select_fun=select_chans_int;
			    break;
			case bin_uint:
			    inp_def.select_fun=select_chans_uint;
			    break;
			case bin_float:
			    inp_def.select_fun=select_chans_float;
			    break;
			case bin_double:
			    inp_def.select_fun=select_chans_double;
			
		    } 	
		    if(inp_def.file_fmt==bin_edf)
			inp_def.no_of_vars=1;

		    inp_def.bytes_per_chan=inp_def.no_of_vars*inp_def.bytes_per_var;
		    inp_def.win_size=inp_def.no_of_chans*inp_def.no_of_vars*inp_def.npoints;
		    inp_def.bytes_win_size=inp_def.nchans*inp_def.bytes_per_chan;
		    inp_def.dwin_size=(inp_def.mtg.mchans+inp_def.mtg.navg)*inp_def.no_of_vars;
		    if(inp_def.mtg.type>MTG_SEL)inp_def.dwin_size+=inp_def.no_of_vars;

	   	 
		    // create input buffers

		    n=inp_def.bytes_win_size*inp_def.npoints;
		    if(inp_def.file_fmt==bin_int12)n+=4; // extra 4 bytes for packed tf

		    if(inp_def.file_fmt==bin_edf)
		    { 
			if((!(inp_def.edf_buf=(short*)malloc(inp_def.bytes_per_win))))
			    mem_error(sname,"inp edf_buf");
			if(!edf_tran &&(!(edf_tran=(short*)malloc(inp_def.bytes_per_win))))
			    mem_error(sname,"inp edf_tran");
			inp_def.edf_sbuf=inp_def.edf_buf;
			inp_def.edf_len=0;
			inp_def.rec_len=inp_def.bytes_per_win/inp_def.bytes_win_size;
		    }
		    if(inp_def.file_fmt!=bin_int12 )
			inp_def.bytes_per_win=inp_def.bytes_win_size;
		    if(inp_def.var_len)
			inp_def.nchvars=inp_def.no_of_chans;
		    else
			inp_def.nchvars=inp_def.no_of_chans*inp_def.no_of_vars;
		    concat(&inp_def);  // concatenates channel and var labels and units 

		    if((!(inp_def.rbuf=malloc(n))))
			mem_error(sname,"inp rbuf");
		    if(inp_def.var_len)
			inp_def.dbuf=NULL; // all montages involving dbuf make no sense with var_len
		    else if(inp_def.mtg.type==NO_MTG && inp_def.file_fmt==bin_double) 
			inp_def.dbuf=(double*)inp_def.rbuf;
		    else 
		    {
			if((!(inp_def.dbuf=
			      (double*)calloc(inp_def.dwin_size*inp_def.npoints,sizeof(double)))))
			    mem_error(sname,"inp dbuf");
		    }

		    memcpy(&inp_def.p_info,&p_info,sizeof(struct pinfo));
 
		    if(no_of_inp<0)
		    {
			struct inpdef *tin;
			if(!(tin=calloc(ni,sizeof(struct inpdef))))
			    mem_error(sname,"new in");
			memcpy(tin,in,(ni-1)*sizeof(struct inpdef));
			if(ni>1)free(in);
			in=tin;
		
		    }
		    else
			if(no_of_inp>1)
			    memcpy(&in[ni-1],&inp_def,sizeof(struct inpdef));
		    parallel_read=check_list(ni-1);
		    if(!parallel_read && (inp_mode==mode_parallel))
		    {
			print_error(WARNING,sname,"not same start time in INPUT sections in parallel read mode",NULL);
			parallel_read=1;
		    }
		    else if(inp_mode==mode_serial)
			parallel_read=0;
		}
	    }
	    else if(s==&data_out_sec)
	    {
		nn=1;
		no++;
		if(no_of_out>0 && no>no_of_out)
		{
		    char erbuf[100];
		    sprintf(erbuf,"%d section %s in ctl %s ignored",no,name,fname);
		    print_error(WARNING,sname,erbuf,NULL);
		    ret=fscanf(inp_ops_fil,"%99[^\n]\n",buf_line);
		    if(ret>0 && copy_ops)fprintf(out_ops_fil,"%s\n",buf_line);
		}
		else
		{ //defaults
		    struct mont *vmtg;
		    memcpy(&out_def,&out_dflt,sizeof(struct outdef));
		    out_def.byte_order=MY_BYTE_ORDER;
		    if(out_def.max_len<=0)
			out_def.max_len=650;
		    if(!out_def.file_fmt)
			out_def.file_fmt=-1; 
		    out_def.path=strdup("./");
		    if(!out_def.templ_dig)
			  out_def.templ_dig=4;
		    out_def.alarm=(FILE*)-1;
		    out_def.all_chans=-1;
		    /*
		      out_def.win_len=-1.;
		      out_def.win_shift=out_def.win_rate=-1.;
		    */
		    out_def.calib=-1; 
		    out_def.start_chan++; //dflt   if not changed it will be 1
		 
		    if(parallel_read)
		    {
			if(out_def.no_of_chans<=0)
			    out_def.no_of_chans=all_chans; // dflt sum of all previous inputs
			out_def.inp_no=last_inp-1;  // last input section for parallel

		    }
		    else
		    {
			last_inp=1;
			out_def.inp_no=0;   // for serial it has to be updated for cur_inp
			if(out_def.no_of_chans<=0)
			    out_def.no_of_chans=in->no_of_chans; // it should be same for all
		    }
		    if(all_chans && no==1)
		    {
			all_chans= out_def.no_of_chans;
			create_all_labels(last_inp);
		    }

		    ret=read_sec(s);
		  
		    if(out_def.all_chans<1)
			out_def.all_chans=all_chans;
		    if(out_def.hdr_name==NULL)
		    {
			print_error(FERROR,sname,"OUTPUT header file name is missing",fname);
		  
		    }

		    if(out_def.file_fmt<0)
		    {
			if(out_def.inp_no<0 || out_def.inp_no>last_inp-1)
			{
			    print_error(WARNING,sname,"OUTPUT wrong inp_no or INPUT section not yet defined ",fname);
			    out_def.inp_no=0;
			 
			}

			out_def.file_fmt=in[out_def.inp_no].file_fmt;
		    }

		    if(out_def.file_fmt==text)
		    {
			print_error(FERROR,sname,"Text format not implemented",fname);
		    }
	  
		    if(out_def.byte_order<0)
		    {
			print_error(FERROR,sname," Unknown OUTPUT byte order",fname);
		 
		    }

		    if(!out_def.title)
		    {
			if(out_def.inp_no>-1 && in[out_def.inp_no].title)
			    out_def.title=strdup(in[out_def.inp_no].title);
			else out_def.title="DATA";
		    }
		  
		    switch(out_def.file_fmt)
		    {
			case bin_char:
			case bin_byte:
			    out_def.bytes_per_var=sizeof(char);
			    break;
			case bin_short:
			case bin_ushort:
			    out_def.bytes_per_var=sizeof(short);
			    break;
			case bin_float: 
			    out_def.bytes_per_var=sizeof(float);
			    break;
			case bin_int: 
			case bin_uint:
			    out_def.bytes_per_var=sizeof(int);
			    break;
			case bin_double: 
			    out_def.bytes_per_var=sizeof(double);
			    break;
			case text:
			    print_error(FERROR,sname,"Text format not implemented",fname);
			case bin_int12:
			    print_error(FERROR,sname,"Telefactor format on output not implemented",fname);
			default:
			    print_error(FERROR,sname,"Not defined or unrecognized file format",fname);
		    }
		    out_def.new_out=out_def.count=0;
		    if(out_def.no_of_vars==0)
		    {
			if(out_def.inp_no>-1 )
			{
			    out_def.var_len=in[out_def.inp_no].var_len;
			    out_def.no_of_vars=in[out_def.inp_no].no_of_vars;
		  
			}
			else
			{
			    out_def.var_len=0;
			    out_def.no_of_vars=1;
			}
		    }
		    else if(out_def.no_of_vars<0)
		    {
			out_def.no_of_vars=-out_def.no_of_vars;
			out_def.var_len=out_def.no_of_vars;
		    }
		    else
			out_def.var_len=0;

		
		    if(out_def.no_of_chans<=0)
		    {
		   
			if(out_def.inp_no<0 || out_def.inp_no>last_inp-1)
			{
			    print_error(FERROR,sname,"OUTPUT wrong Numb_chans",fname);
			  
			}
			out_def.no_of_chans=in[out_def.inp_no].no_of_chans;
		    }
		  
		    if(out_def.win_rate>0.)
		    {
			out_def.win_shift=1.0l/out_def.win_rate;
			out_def.rate_flag=1;
		    }
		    else 
		    {
			if(out_def.inp_no>-1 && out_def.win_shift<=0.)
			    if( in[out_def.inp_no].rate_flag)
			    {
				out_def.win_rate=in[out_def.inp_no].win_rate;
				out_def.rate_flag=1;
				out_def.win_shift=1.0l/out_def.win_rate;
			    }
			    else
			    {
				out_def.win_shift=in[out_def.inp_no].win_shift; // dflt
				out_def.win_rate=1.0l/out_def.win_shift;
				out_def.rate_flag=0;
			    }
			else if(out_def.win_shift>0.)
			{
			    out_def.win_rate=1.0l/out_def.win_shift;
			    out_def.rate_flag=0;
			}  
			else if(out_def.samp_rate>0.)
			{
			    out_def.win_shift=1.0l/out_def.samp_rate;
			    out_def.win_rate=out_def.samp_rate;
			}
			else
			    print_error(FERROR,sname,"OUTPUT Win_shift or Samp_rate has to be defined",fname);
			  
		    }
		  
		    if(out_def.win_len<=0.)
		    { 
			if(out_def.inp_no>-1)
			    out_def.win_len=in[out_def.inp_no].win_len; //dflt
			else
			    out_def.win_len=out_def.win_shift;
		    }  

		    out_def.bytes_per_win=out_def.bytes_per_var*out_def.no_of_vars*out_def.no_of_chans; 
		    if(out_def.inp_no>-1 && out_def.samp_rate<=0)
			out_def.samp_rate=in[out_def.inp_no].samp_rate;
		    if(out_def.inp_no>-1 && out_def.chans_per_file==0)
			if(in[out_def.inp_no].files)
			    out_def.chans_per_file=-1;
			else
			    out_def.chans_per_file=0;
		    if(out_def.chans_per_file==-1)
		    {
			if(!(out_def.files=(FILE**)calloc(out_def.no_of_chans,sizeof(FILE*))))
			    mem_error(sname,"out files");
			if(!(out_def.mbytes=(uint*)calloc(out_def.no_of_chans,sizeof(uint))))
			    mem_error(sname,"out files mbytes");
			if(out_def.var_len)
			    if(!(out_def.no_chan_vars=(int*)calloc(out_def.no_of_chans,sizeof(int))))
				mem_error(sname,"out no_of chans"); 
		    }
		    else
		    {
			out_def.files=NULL;
			out_def.mbytes=NULL;
			out_def.chans_per_file=out_def.no_of_chans;
		    }
		
		    if(!out_def.lst_name)
		    { 
			char *c;
			strcpy(name,out_def.hdr_name);
			c=strchr(name,'.');
			if(c)*c='\0';
			strcat(name,".lst");
			out_def.lst_name=strdup(name);
		    }
		    // template
		    if(out_def.templ)
		    { 
			char *c;
			out_def.ext=strchr(out_def.templ,'.');
			if(out_def.ext!=NULL)
			    *out_def.ext++='\0';
			if((c=strchr(out_def.templ,'#')))
			{
			    out_def.templ_dig=strlen(out_def.templ)-1;
			    out_def.file_no=atoi(c+1);
			    *c='\0';
			    out_def.templ_dig-=strlen(out_def.templ);
			 
			}
		
		    }
		    else
		    { 
			char *c;
			// default templ
			strcpy(name,out_def.hdr_name);
			c=strchr(name,'.');
			if(c)*c='\0';
			out_def.templ=strdup(name);;	  
			out_def.ext=strdup("dat");
		    }
		  
		    if(!out_def.chan_labels)
		    { 
			int ch;
			if(!(out_def.chan_labels=(char**)calloc(out_def.no_of_chans+1,sizeof(char*))))
			    mem_error(sname,"out chan_labels");
			set_lab_no(out_def.chan_labels,out_def.no_of_chans);
			for(ch=0;ch<out_def.no_of_chans;ch++)
			    if(all_labels) out_def.chan_labels[ch+1]=all_labels[ch+out_def.start_chan]; // start_chan is from 1
			    else 
			    {  
				char buf[100];
				sprintf(buf,"Chan_%d",ch+1);
				out_def.chan_labels[ch+1]=strdup(buf);
			    }
		    }
		    if(!out_def.chan_units && all_units)
		    { 
			int ch;
			if(!(out_def.chan_units=(char**)calloc(out_def.no_of_chans+1,sizeof(char*))))
			    mem_error(sname,"out chan_units");
			set_lab_no(out_def.chan_units,out_def.no_of_chans);
			for(ch=0;ch<out_def.no_of_chans;ch++)
			    out_def.chan_units[ch+1]=all_units[ch+out_def.start_chan];
		    }
		    if(out_def.inp_no>-1)
		    {
			vmtg=&in[out_def.inp_no].vmtg;
			if(!out_def.var_labels && vmtg->chan_labels && vmtg->mchans==out_def.no_of_vars)
			{ 
			    int ch;
			    if(!(out_def.var_labels=(char**)calloc(out_def.no_of_vars+1,sizeof(char*))))
				mem_error(sname,"out var_labels");
			    set_lab_no(out_def.var_labels,out_def.no_of_vars);
			    for(ch=1;ch<=out_def.no_of_vars;ch++)
				out_def.var_labels[ch]=vmtg->chan_labels[ch]; 
			}
			if(!out_def.var_units && vmtg->chan_units && vmtg->mchans==out_def.no_of_vars )
			{ 
			    int ch;
			    if(!(out_def.var_units=(char**)calloc(out_def.no_of_vars+1,sizeof(char*))))
				mem_error(sname,"out var_units");
			    set_lab_no(out_def.var_units,out_def.no_of_vars);
			    for(ch=1;ch<=out_def.no_of_vars;ch++)
				out_def.var_units[ch]=vmtg->chan_units[ch];
			}
		    }
		    out_def.start_chan--;   // now it is from 0

		    if(out_def.byte_order != MY_BYTE_ORDER)
			switch(out_def.file_fmt)
			{
			    case bin_short:
			    case bin_ushort:  
				out_def.swap_fun=swap2;
				break;
			    case bin_int:
			    case bin_uint:
			    case bin_float:
				out_def.swap_fun=swap4;
				break;
			    case bin_double:
				out_def.swap_fun=swap8;
			}  
		    out_def.max_len=out_def.max_len*(1<<20);
		  
		    if(out_def.calib<0)
		    {
			if(out_def.inp_no>-1)
			    out_def.calib=1-in[out_def.inp_no].calib; 
			else
			    out_def.calib=0;
		    }
		  

		    switch(out_def.file_fmt)
		    {
			case bin_char:
			    out_def.write_chan_data=write_chan_char_data;
			    out_def.write_data=write_char_data;
			    break;
			case bin_byte:
			    out_def.write_chan_data=write_chan_byte_data;
			    out_def.write_data=write_byte_data;
			    break;
			case bin_short:
			    out_def.write_chan_data=write_chan_short_data;
			    out_def.write_data=write_short_data;
			    break;
			case bin_ushort:  
			    out_def.write_chan_data=write_chan_ushort_data; 
			    out_def.write_data=write_ushort_data;
			    break;
			case bin_int:
			    out_def.write_chan_data=write_chan_int_data;
			    out_def.write_data=write_int_data;	
			    break;
			case bin_uint:
			    out_def.write_chan_data=write_chan_uint_data;
			    out_def.write_data=write_uint_data;
			    break;
			case bin_float:
			    out_def.write_chan_data=write_chan_float_data;
			    out_def.write_data=write_float_data;
			    break;
			case bin_double:
			    out_def.write_chan_data=write_chan_double_data;
			    out_def.write_data=write_double_data;
		    }  

		    if(out_def.path[strlen(out_def.path)-1]!='/')
		    {
			char *tmp=malloc(strlen(out_def.path)+2);
			strcpy(tmp,out_def.path);
			strcat(tmp,"/");
			free(out_def.path);
			out_def.path=tmp;
		    }
		    // actual opening is in first write
		    if(no_of_out<0)
		    {
			struct outdef *tou;
			if(!(tou=calloc(no,sizeof(struct outdef))))
			    mem_error(sname,"new ou");
			memcpy(tou,ou,(no-1)*sizeof(struct outdef));
			if(no>1)free(ou);
			ou=tou;
			memcpy(&ou[no-1],&out_def,sizeof(struct outdef));
		    }
		    else
			if(no_of_out>1)
			    memcpy(&ou[no-1],&out_def,sizeof(struct outdef));

		}
	    }// end data out section
	}// end lookup
	else
	    ret=fscanf(inp_ops_fil,"%99[^\n]\n",buf_line);
    } // end of while loop

    if(no_of_inp>-1 && ni!=no_of_inp)
    {
	print_error(FERROR,sname,"wrong no of  input defs",fname);
  
    }
    else
	no_of_inp=ni;
    if(no_of_out>-1 && no!=no_of_out)
    {
	print_error(FERROR,sname,"wrong no of output defs",fname);
 
    }
    else
	no_of_out=no;
    //rewind(fil);
    del_list[del_no++]="INPUT_OUTPUT";
    del_list[del_no++]="INPUT";
    del_list[del_no++]="OUTPUT";

    if(no_of_inp>0 && !all_labels)
	create_all_labels(no_of_inp);
    if(head && (read_all(fil, head,1)<0))
    {
	print_error(FERROR,sname,"reading tail of ctl file",fname);
 
    };
    fclose(fil);
    return ni;
} 

int check_list(int ni)
{
    // checks if beg in list file and sel file are same for inputs
    // checks type too
    // returns 1 if same as first (parallel reading possible)
    struct inpdef *inp=in+ni;
    char *sname="check_list";
    double beg,begs,end;
    char name[500];
    if(ni==0)
	first_title=strdup(inp->title);
    else
	if(strcmp(inp->title,first_title))
	    return 0;
    if(inp->lst_name && inp->lst_name[0]!='-' && inp->lst_name[0]!='$')
    {
	if(*inp->lst_name=='/'||!strncmp(inp->lst_name,"./",2))
	    strcpy(name,inp->lst_name);
	else
	{
	    strcpy(name,inp->path);
	    strcat(name,inp->lst_name);
	}
  
	inp->lst=fopen(name,"r");
	if(!inp->lst || get_list(inp->lst,inp->freq,&beg,&end))
	{ 
	    print_error(FERROR,sname,"Can't open input list file",name);
	  
	}
//	fclose(inp->lst);inp->lst=NULL;	
    }
    else
	beg=0.;  // pipes are always parallel
    if(inp->sel_name && inp->sel_name[0]!='-' && inp->sel_name[0]!='$')
    {
	if(*inp->sel_name=='/'||!strncmp(inp->sel_name,"./",2))
	    strcpy(name,inp->sel_name);
	else
	{
	    strcpy(name,inp->path);
	    strcat(name,inp->sel_name);
	}  
	inp->sel=fopen(name,"r");
	if(!inp->sel || get_selection(inp->sel,inp->win_shift,&begs,&end,NULL,NULL))
	{ 
	    print_error(FERROR,sname,"Can't open input select file",name);
	  
	}
//	fclose(inp->sel);inp->sel=NULL; 
    }
    else
	begs=beg;
    beg=(begs>beg)?begs:beg;
    if(ni==0)
	beg_first=beg;

    return (beg==beg_first);
 
}

//this is a test

int init_io(int argc,char **argv, struct section *list, int flag)
{
    // initializes all input and output
    // list and select files are open for input and output
    // if not pipe get_epochs is called with global merge_selections flag
    // first input file is open but not output
    // sets f_cur and e_end for first read_next
    // fname is a name of control file
    // list if not NULL points to definitions of control sections
    // flag if not 0  writes out headers based on last inp hdr, mtg and ctl files
    // flag >0 writes out lengths based on INPUT windows params from ctl file
    // flag <0 writes out lengths based on DATA window params from input header

    int i;
    char name[MAX_FILE_NAME],*fname,*sname="init_io";
    struct inpdef *inp;
    struct outdef *out;
 
    progname=argv[0];

#ifdef WINDOWS
	gethostbyname(host);
#else
    gethostname(host,100);
#endif


    if(argc <2 || argc >3 || (argc==3 && strcmp(argv[2],"-r")))
    {
	char *c=strrchr(argv[0],'/');
	if(!c)c=argv[0];
	else c++;
   
	sprintf(name,"%s %s, data_io %s, cfg_io %s\nUsage:%s <ctl file name> [-r]",
                         c,version,dataio_ver,cfg_ver,c);
	puts(name);exit(1);
//print_error(ERROR,sname,name,NULL);
    }

#ifdef linux 
    else if(argc==3)
	stderr=stdout;
#endif

   
    fname=argv[1];
  
    read_cfg(fname,list);  // here header is read twice !

    for(i=0;i<no_of_inp;i++)
    {
	inp=in+i;

	if(inp->lst_name && inp->lst_name[0]!='-')
	{
	    if(*inp->lst_name=='/'||!(strncmp(inp->lst_name,"./",2)))
		strcpy(name,inp->lst_name);
	    else
	    {
		strcpy(name,inp->path);
		strcat(name,inp->lst_name);
	    }
  
	    inp->lst=fopen(name,"r");
	    if(!inp->lst)
	    { 
		print_error(FERROR, sname,"Can't open input list file",name);
	 
	    }
	    inp->fil=NULL;
	    if(inp->lst_name[0]=='$')
		inp->pipe=1; 
	}
	else
	    inp->lst=stdin; //stdin if no list name or name starts with '-'	   

	inp->sel=NULL;
	if(inp->sel_name!=NULL) // optional, if NULL all selected
	{ 
	    if(*inp->sel_name=='/'||!strncmp(inp->sel_name,"./",2))
		strcpy(name,inp->sel_name);
	    else
	    {	
		strcpy(name,inp->path);
		strcat(name,inp->sel_name);
	    }
	    inp->sel=fopen(name,"r");
	    if(!inp->sel)
	    { 
		print_error(FERROR, sname,"Can't open select file",name);
	  
	    }
	}
	else
	{
	    if(inp->epoch)
	    { char time[30],text[200];
	    double len,tod;
	    int id;
 
	    id=0;
	    *text=0;
	    if(sscanf(inp->epoch,"%29[^\t]\t%lf\t%lf\t%d\t%[^\n]\n",time,&len,&tod,&id,text)<3)
	    {		 
		print_error(ERROR,sname,"reading from inp->epoch",NULL);
	    }
 	    tod=rint(tod*inp->freq)/inp->freq;   //rounding to sampling rate
            len=rint(len*inp->freq)/inp->freq;   //rounding to sampling rate
	    inp->e_start=get_time(time,tod);
	    inp->e_end=inp->e_start+len;
	    if(id)
	    {
		inp->sel_id=id;
		inp->sel_text=strdup(text);
		inp->sel_comment=1;
	    }
	    else
		inp->sel_comment=0;
	    }
	    else
	    {
		inp->e_start=0.;	
		inp->e_end=1.e100; //infinity
		inp->sel_comment=0;
		inp->sel_id=0;
		inp->sel_text=strdup("");  // empty string
	    }
	 
	}
	
 	if(!inp->pipe)
	    get_epochs(inp,merge_selections);
	
	if(new_inp_file(inp,1)) // first call 
	{
	    print_error(ERROR, sname,"eof at the begining of  list or select file",NULL);
	}
	if(!inp->pipe)
		{fclose(inp->lst);inp->lst=NULL;}
	if(inp->sel)
	    inp->e_end_pos=-1;// to start
	if(!merge_selections)
		inp->cur_epoch=1;
	inp->gap=1;
	inp->new_out=0;

	

	inp_sel+=(inp->sel!=NULL);
  
    }// end of in loop
  
 
 
 
    for(i=0;i<no_of_out;i++)
    {
	out=ou+i;   
	out->pipe=0;
	if(out->lst_name && out->lst_name[0]!='-')
	{
	    if(*out->lst_name=='/'||!strncmp(out->lst_name,"./",2))
		strcpy(name,out->lst_name);
	    else
	    {
		strcpy(name,out->path);
		strcat(name,out->lst_name);
	    }
#ifndef _WIN32
	    if(out->lst_name[0]=='$')
	    {
		// named pipes
		out->pipe=1;
		unlink(name);
		if(mkfifo(name,S_IRUSR|S_IWUSR))
		{
		    print_error(FERROR, sname,"Can't create fifo lst pipe",name);
		 
		}
		// if lst is pipe files has to be pipe too
	  
	    }
#endif	  
	    if(out->lst_name && out->lst_name[0]!='-')
	    {
		if(*out->lst_name=='/'||!strncmp(out->lst_name,"./",2))
		    strcpy(name,out->lst_name);
		else
		{
		    strcpy(name,out->path);
		    strcat(name,out->lst_name);
		}
		out->fil=out->lst=NULL;
	    }

		   
	}
	else 
	    out->lst=stdout;

	if(*out->hdr_name=='/'||!strncmp(out->hdr_name,"./",2))
	    strcpy(name,out->hdr_name);
	else
	{
	    strcpy(name,out->path);
	    strcat(name,out->hdr_name);
	}

	if(i==0 && out->alarm && out->inp_no>-1 && in[out->inp_no].file_fmt==bin_int12)
	{
	    // alarms from telefactor
	    char *c;
	    c=strrchr(name,'.');
	    if(c)*c='\0';
	    strcat(name,".alarm");
	    out->alarm=fopen(name,"w");
	    if(!out->alarm)
	    { 
		print_error(FERROR, sname,"Can't open output alarm file",name);
	 
	    }
	 
	}
	else
	    out->alarm=NULL;
	out->start=-1.;
	// default for splitting output
/*	if(out->inp_no>-1)
	    out->all_chans=in[out->inp_no].no_of_chans;
	else
	    out->all_chans=out->no_of_chans;
*/
        if(out->join)
	  { int chan;
	    out->join=0;
	  if(out->files)
	    for(chan=0;chan<out->no_of_chans;chan++)
	      new_out_file_chan(out,chan);
	  else
	    new_out_file(out);
	  out->join=1;
	  }
	  
	if(flag)
	{
	    if(flag<0)
		set_out_len(out,in+out->inp_no);
	    write_out_hdr(out);
	    if(list)write_hdr(ctl_name,NULL,list);	  
	    fflush(out->hdr); //live it open for outside use but it has to be closed 
	    // before any i/o operation on other files 
	}

    }// end of ou loop
  
    cur_inp=0;
    inp_eof=0;
    return 0;
}

//Macro defs for get functions, each function
// gets <type> values for one chan numbered like in mtg (from 0) 
// size of data should be n*no_of_vars*sizeof(type)
// n number of points if <1 or >npoints set to npoints
// returns no of elements of <type>


#define GET_DATA_CHAN_DEF(TYPE,CONV)(struct inpdef *inp,TYPE *data,int n,int chan )\
{  \
  int i,j,m,nvars=inp->no_of_vars;  \
  TYPE *d=data;\
  struct mont *mtg=&inp->mtg;  \
  if(n<1)return n; \
  if(n>inp->npoints)n=inp->npoints;\
  if(mtg->type>MTG_SEL)\
 for(i=0,m=0;i<n;i++, m+=inp->dwin_size)   \
 { int l=m+mtg->chan_no[chan]*nvars,k=m+mtg->ref_no[chan]*nvars;   \
   for(j=0;j<nvars;j++)\
 *d++=(TYPE)CONV(inp->dbuf[l++]-inp->dbuf[k++]);   \
 } \
  else if(inp->dbuf)   \
for(i=0,m=0;i<n;i++, m+=inp->dwin_size)\
 { int l=m+mtg->chan_no[chan]*nvars;   \
   for(j=0;j<nvars;j++)\
 *d++=(TYPE)CONV(inp->dbuf[l++]);  \
 } \
  else \
{  char* b=(char*)inp->rbuf+chan*inp->bytes_per_chan; \
if(inp->no_chan_vars) nvars=inp->no_chan_vars[chan];   \
for(i=0;i<n;i++,b+=inp->bytes_per_win,d+=inp->bytes_per_chan)  \
   memcpy(d,inp->rbuf,inp->bytes_per_var*nvars*n); \
}  \
  return n;\
}  
/* 
   int get_double_chan(struct inpdef *inp,double *data,int n,int chan )
   {  
   int i,j,m,nvars=inp->no_of_vars;  
   double *d=data;
   struct mont *mtg=&inp->mtg;  
   if(n<1)return n; 
   if(n>inp->npoints)n=inp->npoints;
   if(mtg->type>MTG_SEL)
   for(i=0,m=0;i<n;i++, m+=inp->dwin_size)   
   { int l=m+mtg->chan_no[chan]*nvars,k=m+mtg->ref_no[chan]*nvars;   
   for(j=0;j<nvars;j++)
   *d++=(double)(inp->dbuf[l++]-inp->dbuf[k++]);   
   } 
   else if(inp->dbuf)   
   for(i=0,m=0;i<n;i++, m+=inp->dwin_size)
   { int l=m+mtg->chan_no[chan]*nvars;   
   for(j=0;j<nvars;j++)
   *d++=(double)(inp->dbuf[l++]);  
   } 
   else 
   {  char* b=(char*)inp->rbuf+chan*inp->bytes_per_chan; 
   if(inp->no_chan_vars) nvars=inp->no_chan_vars[chan];   
   for(i=0;i<n;i++,b+=inp->bytes_per_win,d+=inp->bytes_per_chan)  
   memcpy(d,inp->rbuf,inp->bytes_per_var*nvars*n); 
   }  
   return n;
   }   
*/  
int get_double_chan GET_DATA_CHAN_DEF(double,)
int get_float_chan GET_DATA_CHAN_DEF(float,)
int get_int_chan GET_DATA_CHAN_DEF(int,rint)
int get_uint_chan GET_DATA_CHAN_DEF(uint,rint)
int get_short_chan GET_DATA_CHAN_DEF(short,rint)
int get_ushort_chan GET_DATA_CHAN_DEF(ushort,rint)
int get_char_chan GET_DATA_CHAN_DEF(char,rint)
int get_byte_chan GET_DATA_CHAN_DEF(unsigned char,rint)


   


 


#define GET_DATA_DEF(TYPE,CONV)(struct inpdef *inp,TYPE *data,int n)\
{   \
  int i,j,chan,m,nvars; \
  TYPE *d=data; \
  struct mont *mtg=&inp->mtg;   \
  if(n<1)return n;  \
  if(n>inp->npoints)n=inp->npoints; \
  nvars=inp->no_of_vars;\
  if(mtg->type>MTG_SEL) \
  for(i=0,m=0;i<n;i++, m+=inp->dwin_size)   \
{   \
  for(chan=0;chan<inp->no_of_chans;chan++)  \
	{ int l=m+mtg->chan_no[chan]*nvars,k=m+mtg->ref_no[chan]*nvars;\
	for(j=0;j<nvars;j++)\
	  *d++=(TYPE)CONV(inp->dbuf[l++]-inp->dbuf[k++]);   \
	}   \
}   \
  else if(inp->dbuf)\
for(i=0,m=0;i<n;i++, m+=inp->dwin_size) \
{   \
  for(chan=0;chan<inp->no_of_chans;chan++)  \
	{ int l=m+mtg->chan_no[chan]*nvars; \
	  for(j=0;j<nvars;j++) \
	  *d++=(TYPE)CONV(inp->dbuf[l++]);  \
	}   \
}   \
  else if(inp->no_chan_vars)\
  for(chan=0;chan<inp->no_of_chans;chan++,d+=inp->bytes_per_chan)\
   memcpy(d,inp->rbuf,inp->bytes_per_var*inp->no_chan_vars[chan]);\
  else  \
 memcpy(d,inp->rbuf,inp->bytes_per_win*n);\
  return n; \
}
// get all data assumes parallel reading
#define GET_ALL_DATA_DEF(TYPE,CONV)(TYPE *data,int n)   \
{   \
  int i;\
  TYPE *d=data; \
  if(n<1)return n;  \
  if(n>in->npoints)n=in->npoints;   \
  for(i=0;i<n;i++)  \
{ int ii;   \
  for(ii=0;ii<no_of_inp;ii++)   \
   { struct mont *mtg=&in[ii].mtg;  \
 int m=i*in[ii].dwin_size,chan; \
 if(mtg->type>MTG_SEL)  \
	 for(chan=0;chan<in[ii].no_of_chans;chan++) \
	   { int k,l,j; \
	 l=m+mtg->chan_no[chan]*in[ii].no_of_vars;  \
	 k=m+mtg->ref_no[chan]*in[ii].no_of_vars;   \
	 for(j=0;j<in[ii].no_of_vars;j++)   \
	   *d++=(TYPE)CONV(in[ii].dbuf[l++]-in[ii].dbuf[k++]);  \
	   }\
 else if(in[ii].dbuf)   \
	 for(chan=0;chan<in[ii].no_of_chans;chan++) \
	   { int j,l;   \
	 l=m+mtg->chan_no[chan]*in[ii].no_of_vars;  \
	 for(j=0;j<in[ii].no_of_vars;j++)   \
	   *d++=(TYPE)CONV(in[ii].dbuf[l++]);   \
	   }\
	 else if(in[ii].no_chan_vars)   \
   for(chan=0;chan<in[ii].no_of_chans;chan++,d+=in[ii].bytes_per_chan)\
memcpy(d,in[ii].rbuf,in[ii].bytes_per_var*in[ii].no_chan_vars[chan]);\
 else   \
	   {memcpy(d,in[ii].rbuf,in[ii].bytes_per_win); \
	d=(TYPE*)((char*)d+in[ii].bytes_per_win); } \
   }\
}   \
  return n; \
}

int get_raw_data(struct inpdef *inp,void *data,int n)
{// returns data without conversion(except decoding for telefactor)
   
    int i,j,m; 
    struct mont *mtg=&inp->mtg; 
    char *d=data,*s=(char*)inp->rbuf;
    if(n<1)return n; 
    if(n>inp->npoints)n=inp->npoints;
    if(inp->var_len)n=1;
   
    switch(mtg->type)
    {
	case NO_MTG:
 
	    m=inp->bytes_win_size*n;
	    memcpy(data,inp->rbuf,m);
	    return n;
	case MTG_SEL:
 
	    m=inp->no_of_vars*inp->bytes_per_var;
	    for(j=0;j<n;j++)
	    {
		for(i=0;i<mtg->mchans;i++,d+=m)	 
		    memcpy(d,s+inp->bytes_per_chan*inp->mtg.ch_no[i],m);	   
		s+=inp->bytes_per_win;
	    }  
	    return n;
    }
    print_error(ERROR,"get_raw_data","can't be used with ref montage",NULL);
 

}

int get_raw_all_data(void *data,int n)
{   
    int i,j,ii,m,mm; 
    struct mont *mtg=&in->mtg; 
    char *d=data;
    char *s;
    if(n<1)return n; 
    if(n>in->npoints)n=in->npoints;
    if(in->var_len)n=1;
   
    if(mtg->type>MTG_SEL)
    {
	print_error(ERROR,"get_raw_all_data","can't be used with ref montage",NULL);
  
    }
 
    for(j=0;j<n;j++)
    {
	for(ii=0;ii<no_of_inp;ii++)
	{
	    s=(char*)(in[ii].rbuf)+j*in[ii].bytes_per_win;
	    m=mm=in[ii].no_of_vars*in[ii].bytes_per_var; 
	    if(in[ii].var_len)
	    {
		m=in[ii].var_len*in[ii].bytes_per_var;
	    }
	  
	    for(i=0;i<in[ii].mtg.mchans;i++,d+=mm)	
		if(in[ii].no_chan_vars)
		    memcpy(d,s+in[ii].bytes_per_chan*in[ii].mtg.ch_no[i],in[ii].no_chan_vars[i]*in[ii].bytes_per_var);
		else
		    memcpy(d,s+in[ii].bytes_per_chan*in[ii].mtg.ch_no[i],m);
	}
    }  
    return n;
}


int get_raw_chan(struct inpdef *inp,void *data,int n,int ch)
{
    int j,m;
    struct mont *mtg=&inp->mtg; 
    char *d=data,*s=(char*)inp->rbuf;
    if(n<1)return n; 
    if(n>inp->npoints)n=inp->npoints;
    if(inp->var_len)n=1;
   
    if(mtg->type<MTG_REF)
    {
	s=s+inp->bytes_per_chan*inp->mtg.ch_no[ch];
  
	if(inp->no_chan_vars)
	    m=inp->no_chan_vars[ch]*inp->bytes_per_var;
	else
	    m=inp->no_of_vars*inp->bytes_per_var;
	for(j=0;j<n;j++,d+=m,s+=inp->bytes_per_win)
	    memcpy(d,s,m);
 
	return n;
    }
    print_error(ERROR,"get_raw_chan","can't be used with ref montage",NULL);
 
}


/*
int get_double_dat(struct inpdef *inp,double *data,int n)
{
    int i,j,chan,m,nvars=inp->no_of_vars;
    double *d=data;
    struct mont *mtg=&inp->mtg;
    if(n<1)return n;
    if(n>inp->npoints)n=inp->npoints;

    if(mtg->type>MTG_SEL)
	for(i=0,m=0;i<n;i++, m+=inp->dwin_size)
	{
	    for(chan=0;chan<inp->no_of_chans;chan++)
	    {
		int l=m+mtg->chan_no[chan]*nvars,k=m+mtg->ref_no[chan]*nvars;
		for(j=0;j<nvars;j++)
		    *d++=inp->dbuf[l++]-inp->dbuf[k++];
	    }
	}
    else if(inp->dbuf)
	for(i=0,m=0;i<n;i++, m+=inp->dwin_size)
	{
	    for(chan=0;chan<inp->no_of_chans;chan++)
	    {
		int l=m+mtg->chan_no[chan]*nvars;
		for(j=0;j<nvars;j++)
		    *d++=inp->dbuf[l++];
	    }
	}
    else if(data)
	memcpy(data,inp->rbuf,inp->bytes_per_win*n);
    else
	data=(double*)inp->rbuf;
    return n;

}
*/
/*
  int get_double_all_data(double *data,int n)   
  {   
  int i;
  double *d=data; 
  if(n<1)return n;  
  if(n>in->npoints)n=in->npoints;   
  for(i=0;i<n;i++)  
  { int ii;   
  for(ii=0;ii<no_of_inp;ii++)   
  { struct mont *mtg=&in[ii].mtg;
  int m=i*in[ii].dwin_size,chan;
  if(mtg->type>MTG_SEL)  
  for(chan=0;chan<in[ii].no_of_chans;chan++) 
  { int k,l,j; 
  l=m+mtg->chan_no[chan]*in[ii].no_of_vars;  
  k=m+mtg->ref_no[chan]*in[ii].no_of_vars;   
  for(j=0;j<in[ii].no_of_vars;j++)   
  *d++=(double)(in[ii].dbuf[l++]-in[ii].dbuf[k++]);  
  }
  else if(in[ii].dbuf)   
  for(chan=0;chan<in[ii].no_of_chans;chan++) 
  { int j,l;   
  l=m+mtg->chan_no[chan]*in[ii].no_of_vars;  
  for(j=0;j<in[ii].no_of_vars;j++)   
  *d++=(double)(in[ii].dbuf[l++]);   
  }
  else if(in[ii].no_chan_vars)   
  for(chan=0;chan<in[ii].no_of_chans;chan++,d+=in[ii].bytes_per_chan)
  memcpy(d,in[ii].rbuf,in[ii].bytes_per_var*in[ii].no_chan_vars[chan]);
  else   
  {memcpy(d,in[ii].rbuf,in[ii].bytes_per_win); 
  d=((char*)d+in[ii].bytes_per_win); }   
  }
  }   
  return n; 
  }
*/
int get_double_dat GET_DATA_DEF(double,)
int get_float_data GET_DATA_DEF(float,)
int get_int_data GET_DATA_DEF(int,rint)
int get_uint_data GET_DATA_DEF(uint,rint)
int get_short_data GET_DATA_DEF(short,rint)
int get_ushort_data GET_DATA_DEF(ushort,rint)
int get_char_data GET_DATA_DEF(char,rint)
int get_byte_data GET_DATA_DEF(unsigned char,rint)

int get_double_all_data GET_ALL_DATA_DEF(double,)
int get_float_all_data GET_ALL_DATA_DEF(float,)
int get_int_all_data GET_ALL_DATA_DEF(int,rint)
int get_uint_all_data GET_ALL_DATA_DEF(uint,rint)
int get_short_all_data GET_ALL_DATA_DEF(short,rint)
int get_ushort_all_data GET_ALL_DATA_DEF(ushort,rint)
int get_char_all_data GET_ALL_DATA_DEF(char,rint)
int get_byte_all_data GET_ALL_DATA_DEF(unsigned char,rint)

int get_double_data_mtg(struct inpdef *inp, double *data, int n, struct mont *mtg)
{
// allows for remontage of data
// original montage is lost
    inp->mtg=*mtg;
    return get_double_dat(inp,data,n);
}

int get_double_data(struct inpdef *inp, double *data, int n)
{
    if((inp->mtg).type<=MTG_SEL && inp->dbuf)
    {
	if(data)
	    memcpy(data,inp->dbuf,n*inp->dwin_size*sizeof(double));
	else
	    data=inp->dbuf;
	return n;
    }
    else
	return get_double_dat(inp,data,n);
}

//reads and gets data, returns number get, if it is less than n returns negative value or 0



int read_double_data(struct inpdef *inp,double *data,int n,int part) 
{
    return get_double_data(inp,data,read_next_win(inp,n,part));
}

int read_float_data(struct inpdef *inp,float *data,int n,int part) 
{
    return get_float_data(inp,data,read_next_win(inp,n,part));
}

int read_int_data(struct inpdef *inp,int *data,int n,int part) 
{
    return get_int_data(inp,data,read_next_win(inp,n,part));
}

int read_uint_data(struct inpdef *inp,uint *data,int n,int part) 
{
    return get_uint_data(inp,data,read_next_win(inp,n,part));
}

int read_short_data(struct inpdef *inp,short *data,int n,int part) 
{
    return get_short_data(inp,data,read_next_win(inp,n,part));
}

int read_ushort_data(struct inpdef *inp,ushort *data,int n,int part) 
{
    return get_ushort_data(inp,data,read_next_win(inp,n,part));
}

int read_char_data(struct inpdef *inp,char *data,int n,int part) 
{
    return get_char_data(inp,data,read_next_win(inp,n,part));
}

int read_byte_data(struct inpdef *inp,unsigned char *data,int n,int part) 
{   
    return get_byte_data(inp,data,read_next_win(inp,n,part));
}

// reads and gets data from all input streams

 
int read_raw_all_data(void *data,int n,int part) 
{
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_raw_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof)
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_raw_data(in+cur_inp,data,nn);
    }
}

int read_double_all_data(double *data,int n,int part)
{ 
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    { 
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_double_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof)
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_double_data(in+cur_inp,data,nn);
    }
}

int read_int_all_data(int *data,int n,int part)
{ 
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    { 
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_int_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof)
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_int_data(in+cur_inp,data,nn);
    }
}

int read_uint_all_data(uint *data,int n,int part)
{
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_uint_all_data(data,nn);
    } 
    else
    { 
	if(in[cur_inp].eof) 
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;   
	return get_uint_data(in+cur_inp,data,nn);
    } 
}
int read_short_all_data(short *data,int n,int part)
{ 
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    { 
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_short_all_data(data,nn);
    } 
    else
    { 
	if(in[cur_inp].eof)
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_short_data(in+cur_inp,data,nn);
    }
}
int read_ushort_all_data(ushort *data,int n,int part)
{
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_ushort_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof)
	    cur_inp++; 
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_ushort_data(in+cur_inp,data,nn);
    }
}

int read_char_all_data(char *data,int n,int part) 
{
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_char_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof) 
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;   
	return get_char_data(in+cur_inp,data,nn);
    }
}


int read_byte_all_data(unsigned char *data,int n,int part)
{
    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    nn=read_next_win(in+i,n,part);
	    inp_eof+=in[i].eof;
	}
	return get_byte_all_data(data,nn);
    }
    else
    {
	if(in[cur_inp].eof)
	    cur_inp++;
	nn=read_next_win(in+cur_inp,n,part);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return get_byte_data(in+cur_inp,data,nn);
    }
}


// Warning: data is converted in place 
#define WRITE_DATA_DEF(TYPE,CONV)(struct outdef *out,double *data,int npoints)  \
{  \
 TYPE* f;  \
 double *d,*dd,*beg;   \
 int i,j,nd,nvars; \
 if(npoints<1)return;  \
 nvars=out->no_of_vars;\
 nd=nvars*out->all_chans;  \
 beg=dd=data+out->start_chan*nvars;\
 f=(TYPE*)beg; \
 for(j=0;j<npoints;j++,dd+=nd) \
 for(i=0,d=dd;i<nvars*out->no_of_chans;i++)\
	 *f++=(TYPE)CONV(*d++);\
 write_data_raw(out,beg,npoints,nvars,sizeof(TYPE));  \
}

/*

  void write_short_data(struct outdef *out,double *data,int npoints)
  {  
  short* f;  
  double *d,*dd,*beg;
  int i,j,nd,nvars;
  if(npoints<1)return;  
  nvars=out->no_of_vars; 
  nd=nvars*out->all_chans;  
  beg=dd=data+out->start_chan*nvars;
  f=(short*)beg;
  for(j=0;j<npoints;j++,dd+=nd)   
  for(i=0,d=dd;i<nvars*out->no_of_chans;i++)   
  *f++=(short)rint(*d++);
  write_data_raw(out,data,npoints,nvars,sizeof(short));   
  }
*/

void write_float_data WRITE_DATA_DEF(float,)
void write_int_data WRITE_DATA_DEF(int,rint)
void write_uint_data WRITE_DATA_DEF(uint,rint)	 
void write_short_data WRITE_DATA_DEF(short,rint)
void write_ushort_data WRITE_DATA_DEF(ushort,rint)
void write_char_data WRITE_DATA_DEF(char,rint)
void write_byte_data WRITE_DATA_DEF(unsigned char,rint)

void write_double_data(struct outdef *out,double *data,int npoints)
{ 
    write_data_raw(out,data,npoints,out->no_of_vars,sizeof(double));
}


void write_all_data_raw(void *data,int npoints,int nvars,int varsize)
{
    int i;   
    if(npoints<1)return; 
    for(i=0;i<no_of_out;i++)
    { 
	write_data_raw(ou+i,data,npoints,nvars,varsize);
    }
}
void write_all_data(double *data,int npoints)
{
    int i;   
    if(npoints<1)return; 
    for(i=0;i<no_of_out;i++)
	ou[i].write_data(ou+i,data,npoints);
}


  

#define WRITE_CHAN_DATA_DEF(TYPE,CONV)(struct outdef *out,double *data,int n,int chan,int nv) \
{ \
 TYPE* f=(TYPE*)data; \
 double *d=data,*dn;  \
 if(n<1)return;   \
 if(nv<0)nv=n=out->no_of_vars;\
 dn=data+n;   \
 while(d<dn)  \
	 *f++=(TYPE)CONV(*d++);   \
 write_chan_data(out,(TYPE*)data,n,chan,nv);  \
}


void write_chan_double_data(struct outdef *out,double *data,int n,int chan,int nv)
{ write_chan_data(out,data,n,chan,nv);}

void write_chan_float_data WRITE_CHAN_DATA_DEF(float,)
void write_chan_int_data WRITE_CHAN_DATA_DEF(int,rint)
void write_chan_uint_data WRITE_CHAN_DATA_DEF(uint,rint)	 
void write_chan_short_data WRITE_CHAN_DATA_DEF(short,rint)
void write_chan_ushort_data WRITE_CHAN_DATA_DEF(ushort,rint)
void write_chan_char_data WRITE_CHAN_DATA_DEF(char,rint)
void write_chan_byte_data WRITE_CHAN_DATA_DEF(unsigned char,rint)

void new_inp_size(struct inpdef *inp, int n)
{
    int flag=(inp->dbuf && inp->dbuf != (double*)inp->rbuf);
    char *sname="new_inp_size";
    if(flag)free(inp->dbuf);
    if(inp->rbuf)free(inp->rbuf);
    if(inp->npoints==inp->spoints)
    {   char err[50];
    inp->spoints=n;
    sprintf(err,"npoints changed  to %d points",n);
    print_error(WARNING,sname,err,NULL);
    inp->win_shift=inp->spoints/inp->freq;
    }
    inp->npoints=n;
    inp->win_len=inp->npoints/inp->freq;
    n=inp->bytes_win_size*inp->npoints;
    if(inp->file_fmt==bin_int12)n+=4; // extra 4 bytes for packed tf
    if((!(inp->rbuf=malloc(n))))
	mem_error(sname,"new inp rbuf");
    if(!flag)inp->dbuf=inp->rbuf; 
    else if((!(inp->dbuf=(double*)calloc(inp->dwin_size*inp->npoints,sizeof(double)))))
	mem_error(sname,"new inp dbuf");
    return;
}

void set_out_fmt(struct outdef *out, int nchans, int nvars, int fmt, int files)
{

    // allows for changes in out parameters
    // 
    char* sname= "set_out_fmt";

  
    if(fmt>=0)
	out->file_fmt=fmt;
  
    if(nvars) 
	out->no_of_vars=nvars;
  
    if(nchans>0)
	out->no_of_chans=nchans;

    if(files<0)
	out->chans_per_file=-1;
    else if(files>0)
	out->chans_per_file=out->no_of_chans;


    switch(out->file_fmt)   
    {
	case bin_char:
	case bin_byte:
	    out->bytes_per_var=sizeof(char);
	    break;
	case bin_short:
	case bin_ushort:
	    out->bytes_per_var=sizeof(short);
	    break;
	case bin_float: 
	    out->bytes_per_var=sizeof(float);
	    break;
	case bin_int: 
	case bin_uint:
	    out->bytes_per_var=sizeof(int);
	    break;
	case bin_double: 
	    out->bytes_per_var=sizeof(double);
	    break;
	case text:
	    print_error(ERROR,sname,"Text format not implemented",NULL);
	case bin_int12:
	    print_error(ERROR,sname,"Telefactor format on output not implemented",NULL);
	default:
	    print_error(ERROR,sname,"Not defined or unrecognized file format",NULL);

    }
 
    if(out->no_of_vars<0)
    {
	out->no_of_vars=-out->no_of_vars;
	out->var_len=out->no_of_vars;
    }
    else
	out->var_len=0;

    out->bytes_per_win=out->bytes_per_var*out->no_of_vars*out->no_of_chans;
}
 
void concat(struct inpdef *inp)
{
    // concatenates labels and units from chans and vars
    int i,nc,nv,j,k,kk;
    char **l[2],**lc[2],**lv[2];
    char c[2]={'_',' '};
    char buf[MAX_LINE];

    nc=inp->no_of_chans;

    if(inp->var_len)
	nv=1;
    else
	nv=inp->no_of_vars; 
    if(nc*nv>1000)
    {
	inp->labels=inp->mtg.chan_labels;
	inp->units=inp->mtg.chan_units;
	return;
    }
    lc[0]=inp->mtg.chan_labels;
    lc[1]=inp->mtg.chan_units;
    lv[0]=inp->vmtg.chan_labels;
    lv[1]=inp->vmtg.chan_units;

    if(nv==1 && lv[0] && !strncmp(lv[0][1],"V1",2)) 
    {
	l[0]=lc[0];
	if(lv[1])l[1]=lv[1];
	else l[1]=lc[1];
    }
    else if(nc==1 && lc[0] && !strncmp(lc[0][1],"C1",2)) 
    {   
	l[0]=lv[0];
	if(lc[1])l[1]=lc[1];
	else l[1]=lv[1];
    }
    else

	for(i=0;i<2;i++)
	{

	    if(!lc[i] && !lv[i])
	    {l[i]=NULL;continue;}

	    l[i]=(char**)calloc((nc*nv+1),sizeof(char*));
	    if(!l[i])
		mem_error("concat","l");
	    set_lab_no(l[i],nc*nv);

	    if(lc[i] && lv[i])
	    {
		for(j=1,kk=1;j<=nc;j++)
		    for(k=1;k<=nv;k++)
		    {
			sprintf(buf,"%s%c%s",lc[i][j],c[i],lv[i][k]);
			l[i][kk++]=strdup(buf);
		    }
	    }
	   
	    else if (!lc[i] && lv[i])
	    {
		for(j=1,kk=1;j<=nc;j++)
		    for(k=1;k<=nv;k++)
			l[i][kk++]=lv[i][k];
	    }
	    else //(lc[i] && ! lv[i])
	    {
		  
		for(j=1,kk=1;j<=nc;j++)
		    for(k=1;k<=nv;k++)
			l[i][kk++]=lc[i][j];
	    }
	}
    inp->labels=l[0];
    inp->units=l[1];

}

int get_list(FILE* fil,double f,double *beg, double *end)
{
  
    double len,tod;
    char time[30];
	
    do{
	if(fscanf(fil,"%*s\t%29[^\t]\t%lf\t%lf%*[^\n]\n",time,&len,&tod)!=3)
	    return 1;
    }while(len<=0);  // skiping lines with negative lengths
    tod=rint(tod*f)/f;   //rounding to sampling rate
    len=rint(len*f)/f;   //rounding to sampling rate
    *beg=get_time(time,tod);
    *end=*beg+len;
    return 0;
}

int get_selection(FILE* fil,double f, double *beg, double *end, int* type_no, char **type)
{
    // reads line from file fil
    // if type_no==NULL skips negative or zero length selections
    // if type not NULL reads comment to type
    // dt is samp interval for rounding
 
  
    double len,tod;
    int no,skip=(type_no==NULL);
    char time[30];
    char buf[MAX_LINE];

 
    if(!type)
	type_no=NULL;

    do{
	*buf=no=0;
	if(fscanf(fil,"%29[^\t]\t%lf\t%lf\t%d\t%[^\n]\n",time,&len,&tod,&no,buf)<3)
	    return 1;
   
    }while((skip && len<=0) || (time[0]=='#'));  // skiping events with negative or zero lengths or comments
    tod=rint(tod*f)/f;   //rounding to sampling rate
    len=rint(len*f)/f;   //rounding to sampling rate
    *beg=get_time(time,tod);
    *end=*beg+len;

    if(type)
    {
	if(*buf==0)
	
	    *type=NULL;
	else
   
	    *type=strdup(buf); 
    }
   
    if(type_no)*type_no=no;
   
    return 0;
}


double *get_epochs(struct inpdef *inp, int merge_sel)
{
    // returns list of intervals common for list end select
    // if merge_sel == 0 select is not used (i.e. only list of epochs in list file)
    // it returns beg and end of epochs
    // in fnames returns names of files
    // format for list is same as in intervals_in
    // format for fnames is same as in labels_in
    // memory is allocated in this procedure
    // it closes both list and select files !!
 
    int n=0,k,i,j,m;
    double *sel,*start;
    char** selname,**sel_text,**tmp_sel_text,*sname="get_epochs";
    char time[30],name[200];
    double len,tod;
    double beg1,end1;
    double *beg,*end,bg,ed,old_beg=0,old_end=0;
    int id,*sel_id,nsel;
    char text[200];
    char *txt;
    double eps=1./inp->samp_rate+1.e-6;  // 1.e-6 is insignificant for time scale with full secs
    if(!inp->lst)goto err_lst;
    if(inp->pipe)
	print_error(FERROR,"get_epochs","this mode not available for pipes",inp->sel_name);
    // count lines in list file
    rewind(inp->lst);
    while(fgets(buf_line,MAX_LINE,inp->lst))
	{
	        if(buf_line[0]=='#' && !inp->zeros) continue;
		n++;  
	}

    if(!feof(inp->lst)) goto err_lst;
   
    m=2*n+1;
    rewind(inp->lst);
    if(!(inp->list=(double*)malloc(m*sizeof(double))))mem_error(sname,"list");
    if(!(inp->start_file=(double*)malloc((n+1)*sizeof(double))))mem_error(sname,"start_file");
    if(!(inp->fname=(char**)malloc((n+1)*sizeof(char*))))mem_error(sname,"fname");
    if(!(inp->id=(int*)malloc((n+1)*sizeof(int))))mem_error(sname,"id");
    if(!(inp->text=(char**)malloc((n+1)*sizeof(char*))))mem_error(sname,"text");

    inp->list[0]=inp->start_file[0]=inp->id[0]=n;
    set_lab_no(inp->fname,n);
    set_lab_no(inp->text,n);
    for(i=0,j=1;i<n;i++,j+=2)
    { 
	
        do{ 
        *text=0;id=0;
	if((!fgets(buf_line,MAX_LINE,inp->lst) ||
	    (sscanf(buf_line,"%s\t%29[^\t]\t%lf\t%lf\t%d\t%[^\n]\n",name,time,&len,&tod,&id,text)<4)) &&
	   !feof(inp->lst)) goto err_lst;}
 	while(buf_line[0]=='#' && !inp->zeros);
        tod=rint(tod*inp->freq)/inp->freq;   //rounding to sampling rate
        len=rint(len*inp->freq)/inp->freq;   //rounding to sampling rate
	inp->list[j]=get_time(time,tod);
	if(fabs(inp->list[j]-inp->list[j-1])<=eps) // equal due to rounding errors in Stellate LOG file
	    inp->list[j]=inp->list[j-1];  // to prevent artifficial gaps between files
	inp->list[j+1]=inp->list[j]+len;
	if(!(inp->fname[i+1]=strdup(name)))mem_error(sname,"fname[i]");
	inp->start_file[i+1]=inp->list[j];
	inp->id[i+1]=id;
        if(strlen(text)<1)
           inp->text[i+1]=NULL;
        else 
          if(!(inp->text[i+1]=strdup(text)))mem_error(sname,"text[i]");
  
  
    }
  
    if(!inp->sel || !merge_sel)
    {
//	fclose(inp->lst);
//	inp->lst=NULL;
	rewind(inp->lst);
	return inp->list;
    }

    // combine selections with list


 
//    first read and combine overlapping selections
    rewind(inp->sel);k=0;
    while(fgets(buf_line,MAX_LINE,inp->sel))k++;
    if(!feof(inp->sel)) goto err_sel;
    DMEM(beg,k);
    DMEM(end,k);
    if(!(tmp_sel_text=(char**)calloc(k,sizeof(char*))))mem_error(sname,"tmp_sel_text");
    merge_selections=0;  //temp change to read selections
    rewind(inp->sel);k=0;
    get_selection(inp->sel,inp->win_shift,&old_beg,&old_end,&id,&txt);
    beg[0]=old_beg;
    if(txt && strlen(txt)>0)
	tmp_sel_text[0]=strdup(txt);
    while(!get_selection(inp->sel,inp->win_shift,&bg,&ed,&id,&txt))
    {
	if(bg<old_beg)
	  print_error(FERROR,sname,"selections_out_of_order",inp->sel_name);
	if(bg>old_end)
          {
		
		end[k++]=old_end;
		beg[k]=old_beg=bg;
		 if(txt && strlen(txt)>0)
			tmp_sel_text[k]=strdup(txt);
          }
	if(ed>old_end)old_end=ed;
    }
    end[k++]=old_end;
    nsel=k;
    j=2*(k+n)+1;  // this is maximum possible from combination even when sel is empty (k==0)
    if(!(sel=(double*)malloc(j*sizeof(double))))mem_error(sname,"sel");
    if(!(start=(double*)malloc((k+n+1)*sizeof(double))))mem_error(sname,"start");
    if(!(selname=(char**)malloc((k+n+1)*sizeof(char*))))mem_error(sname,"selname");
    if(!(sel_id=(int*)malloc((k+n+1)*sizeof(int))))mem_error(sname,"id");
    if(!(sel_text=(char**)calloc((k+n+1),sizeof(char*))))mem_error(sname,"sel_text");
    
    k=j;
    i=j=1;
    id=1;
    beg1=beg[0];
    end1=end[0];
    for(;;) 
    {
	while(i<m && beg1>=inp->list[i])
	{
	    i+=2;
	}
	if(i>m)goto end;
	if(i>1 && beg1<inp->list[i-1])
	{
	    sel[j++]=beg1;
	    selname[j/2]=strdup(inp->fname[i/2]);
	    start[j/2]=inp->start_file[i/2];
            sel_id[j/2]=inp->id[i/2];
	    if(n>k && inp->text[i/2])sel_text[j/2]=strdup(inp->text[i/2]);
            else  sel_text[j/2]=tmp_sel_text[id-1];
	    if(end1>inp->list[i-1])
	    {
		sel[j++]=beg1=inp->list[i-1];
/*		beg1=inp->list[i];
		i+=2;	 
*/
		if(i==m)goto end;
		else continue;

	    }
	    else
		sel[j++]=end1;
	}
	else if(end1>inp->list[i])
	{
	    sel[j++]=inp->list[i];
	    selname[j/2]=strdup(inp->fname[(i+1)/2]);
	    start[j/2]=inp->start_file[(i+1)/2];
	    sel_id[j/2]=inp->id[(i+1)/2];
	    if(inp->text[(i+1)/2]) sel_text[j/2]=strdup(inp->text[(i+1)/2]);
            else sel_text[j/2]=tmp_sel_text[id-1];

	    if(end1>inp->list[i+1])
	    {
		sel[j++]=inp->list[i+1];
		i+=2;
		if(i>m)goto end;
		else 
		{
		    beg1=inp->list[i];	
		    continue;
		}
	    }
	    else
		sel[j++]=end1;
	}
        if(id==nsel) goto end;
        beg1=beg[id];end1=end[id++];
        
    }
 end:
    
    free(inp->list);
    free(inp->start_file);
    for(i=1;i<=n;i++)
    {
	free(inp->fname[i]);
	if(inp->text[i])free(inp->text[i]);
    }
    free(inp->fname);
    free(inp->id);
    free(inp->text);
    rewind(inp->sel);
    sel[0]=start[0]=sel_id[0]=j/2;
    set_lab_no(selname,j/2);
    set_lab_no(sel_text,j/2);
     
    if(j<k) 
    {
	// save some memory
	n=j/2+1;
	if(!(inp->list=(double*)malloc(j*sizeof(double))))mem_error(sname,"2nd list");
	if(!(inp->fname=(char**)malloc(n*sizeof(char*))))mem_error(sname,"2nd fname");
        
	if(!(inp->start_file=(double*)malloc(n*sizeof(double))))mem_error(sname,"2nd start_file");
        if(!(inp->id=(int*)malloc(n*sizeof(int))))mem_error(sname,"2nd id");
        if(!(inp->text=(char**)malloc(n*sizeof(char*))))mem_error(sname,"2nd text");
   
	memcpy(inp->list,sel,j*sizeof(double));  
	free(sel);
	memcpy(inp->start_file,start,n*sizeof(double));  
	free(start);
	memcpy(inp->fname,selname,n*sizeof(char*));
	free(selname);
        memcpy(inp->id,sel_id,n*sizeof(int));
	free(sel_id);
        memcpy(inp->text,sel_text,n*sizeof(char*));
	free(sel_text);
        
    }
    else
    {
	inp->list=sel;
	inp->fname=selname;
        inp->id=sel_id;
        inp->text=sel_text;
    }
    fclose(inp->sel);
    inp->sel=NULL;
    merge_selections=merge_sel;
    inp->cur_epoch=0;
//    inp->e_start=inp->list[1];
//    inp->e_end=inp->list[2]; 
//    fclose(inp->lst);
//    inp->lst=NULL;
   
      rewind(inp->lst);
    free(beg);
    free(end);
    free(tmp_sel_text);
    return inp->list;

 err_sel:
    print_error(FERROR,"get_epochs","reading from  select file",inp->sel_name);
 

  
 err_lst:
    print_error(FERROR,"get_epochs","reading from list file",inp->lst_name);
  
}



int pos_all_files(int j,double sec)
{ 
    int i,ret;
    if(!parallel_read) return pos_file(in+cur_inp,j,sec);
    ret=pos_file(in,j,sec);
    for(i=1;i<no_of_inp;i++)
	if(pos_file(in+i,j,sec)!=ret)
	    print_error(ERROR,"pos_all_files","Not all files in same position",NULL);
    return ret;
}


int pos_file(struct inpdef *inp,int j, double sec)
{
    // positions file at sec in j-th epoch(starting at 1
    // sets inp->cur_epoch;
    // if sec outside of file returns 0
    // otherwise returns j=cur_epoch 
  
    double* list=inp->list+1,last=inp->f_start;
    int i;
    if(j>inp->list[0]||j<1)return 0;
    i=2*j-2;
    if(strncmp(inp->name,inp->fname[j],NAME_LEN))
    {
	strncpy(inp->name,inp->fname[j],NAME_LEN);
	inp->f_start=inp->start_file[j]; 
	inp->f_end=list[i+1];
	if(!inp->sel_comment)
	{
	    inp->sel_id=inp->id[j];
	    inp->sel_text=inp->text[j];
	}
	new_inp_file(inp,0); // no read of lst file	 
    }
    if(merge_selections)
    {
	// definitions of epochs differ from next_win
	// if merged selections
	inp->e_start=list[i];
	inp->e_end=list[i+1];
	if(sec<inp->e_start || sec >=inp->e_end)
	    return inp->cur_epoch=0;
        if(!inp->sel_comment)
	{
	    inp->sel_id=inp->id[i/2+1];
	    inp->sel_text=inp->text[i/2+1];
	}
    }
    else
       if(sec<inp->f_start || sec >=inp->f_end)
	            return inp->cur_epoch=0;
    if(inp->e_start>0.) 
      inp->e_start_pos=irint((inp->e_start-inp->f_start)*inp->freq);
    
    if(inp->e_end<1e100)
      inp->e_end_pos=inp->e_start_pos+irint((inp->e_end-inp->e_start)*inp->freq); 
    else
      inp->e_end_pos=INT_MAX;

    inp->cur+=irint((last-inp->f_start)*inp->freq);
    i=irint((sec-inp->f_start)*inp->freq);
    inp->skip_fun(inp,i);
    inp->cur_epoch=j; // starts at 1
    return j;
}


int find_in_list_all(double sec,double len)
{ // it needs to find it all to position files
  // could be improved by skiping search
    int i,ret;
    ret=find_in_list(in+cur_inp,sec,len);
    if(!parallel_read)return ret;
    for(i=1;i<no_of_inp;i++)
	if(find_in_list(in+i,sec,len)!=ret)
	    print_error(ERROR,"find_in_list_all","Not all streams found the same epoch. \
May be caused by different selection files.",NULL);
    return ret;
}
  

int find_in_list(struct inpdef *inp,double sec,double len)
{
    // finds sec in list and returns index (starting at 1)
    // if in gap returns negative index to previous  
    // if before first returns 0
    // do not reposition list and select files
    // positions stream at sec, or len before end of previous (or begining of first)
    // if len ==0 always positions at begining of next unless after last when it is undef
    // uses info in list,fname,start_file as computed  by get_epochs
    // if continuous epochs k is actual epoch in which sec-len is located
    // relevant only for gaps

    int i,j,k,n=inp->list[0];
    double* list=inp->list+1;

    // search for sec in list
    // assumes that starting time in list is ordered
    // and ending time is ordered but may overlap
    // starts search from cur_epoch forward or backward
 
    j=2*inp->cur_epoch-2;
    // first check the most common 
 
    if(sec<list[0])
	k=0; // before first
    else if(sec>=list[2*n-1])
	k=-n; // after last
    else
    {
	int l=0,h=2*n-1;
	while(h-l>1)
	{
	    i=(l+h)/2;
	    if(sec<list[i])h=i;
	    else l=i;
	}
// list[l]<=sec<list[r] 
	if(l%2==0)
	    k=l/2+1; // in k-th epoch
	else
	    k=-h/2;   // in gap after k-th epoch
    }
	  

    /*  
	while(h-l>1)
	{
	i=(l+h)/2;
	if(sec<list[i])h=i;
	else l=i;
	}
// list[l]<=sec<list[h] 
if(l%2==0)
k=l/2+1; // in k-th epoch
else
k=-h/2;   // in gap after k-th epoch
}
    */
 found:
    i=2*abs(k)-2;
    j=abs(k);	  
    if (k<=0)
    {
	if(len>0)
	{
	    if(k==0)
	    {sec=list[0];j=1;}
	    else
	    {
		sec=list[i+1]-abs(len);
		if(sec<list[i])sec=list[i];
	    }
	}
	else if(j<n)
	{
	    sec=list[i+2];
	    if(sec >= inp->e_end && !merge_selections)
	    {
		// next selection epoch
		if(read_select_line(inp))return -n; // eof
		else 
		{
		    inp->gap=1;
		    return find_in_list(inp,inp->e_start,0);
		}
	    }
	    j++;
	}
	else
	    return k;
    }
    // if in epoch sec do not change   
 
    if(pos_file(inp,j,sec))
    {
	//	   printf("f_start=%12lf,f_end=%12lf,e_start=%12lf,e_end=%12lf,sec=%12lf,j=%d\n",
	//		  inp->f_start,inp->f_end,inp->e_start,inp->e_end,sec,j);
	  
	return k;
    }
    else 
    {
	//	 printf("f_start=%lf,f_end=%12lf,e_start=%12lf,e_end=%12lf,sec=%12lf,j=%d\n",
	//		  inp->f_start,inp->f_end,inp->e_start,inp->e_end,sec,j);
	print_error(FERROR,"find_in_list","pos_file returns 0",inp->sel_name);
	 
    }
	   
}


int read_all_current(int n)
{

    int i,nn;
    inp_eof=0;
    if(parallel_read)
    {
	for(i=0;i<no_of_inp;i++)
	{
	    n=read_current(in+i,n);
	    inp_eof+=in[i].eof;
	}
	return n;
    }
    else
    {
	nn=read_current(in+cur_inp,n);
	if(in[cur_inp].eof && cur_inp==no_of_inp-1)
	    inp_eof=1;
	return nn;
    }
}

int read_current(struct inpdef *inp,int npoints)
{
    // reads np points(windows) starting at pos
    // as set previously by find_in_list
    // does not check if in epoch
    // selects channels into dbuf
    // works for fixed length (for var_len & npoints=1 it also should work)
    // returns no of windows(points) read
    // updates f_cur = begining of current window
    // updates pos,end  and cur  in no of windows
    // always reads partial data
    // if eof on stream inp->eof=1; and 0 read
    // if gap during read sets gap and returns partial
    // new_out and new_epoch are not changed
    // on input inp->cur should be the last position
    // and inp->end should be the end pos of last read
 
    int np,n,k=0;
    char *sbuf;
    double *dbuf;
    int flag=(inp->dbuf && inp->dbuf != (double*)inp->rbuf);

    if(npoints<=0 || npoints>inp->npoints)npoints=inp->npoints;
  
    np=npoints;
    inp->eof=inp->gap=0;
    sbuf=inp->rbuf;
    dbuf=inp->dbuf;
 
    inp->cur=inp->pos;
    inp->f_cur=inp->f_start+inp->cur/inp->freq;
  
    do
    { 
	sync_cur=inp->pos; // this is needed for int12
	n=0;
	if(inp->pos<inp->e_end_pos && (n=inp->read_fun(inp,sbuf,np))>0)
	{
	    if(inp->pos+n>=inp->e_end_pos)
		n=inp->e_end_pos-inp->pos;
	    inp->pos+=n;
	    k+=n;
	    if(inp->swap_fun!=NULL)
		inp->swap_fun((short*)sbuf,inp->bytes_win_size,n);
	    if(flag)
	    {
		inp->select_fun(inp,sbuf,dbuf,n);
		dbuf+=n*inp->dwin_size;
	    }
	    sbuf+=n*inp->bytes_win_size;	 
	}
	else
	{
	 
	    int nn=inp->cur_epoch;

	    if(nn++==inp->list[0])
	    {
		inp->eof=1;
		inp->end=inp->cur+k;
		return k;
	    }
	    else if(inp->e_end!=inp->list[2*nn-1])
	    {  
		inp->gap=1;
		inp->end=inp->cur+k;  
		return k;
	    }
	    pos_file(inp,nn,inp->e_end);// new file but continous	
	}	  
	np-=n; 
	//  pos=inp->pos;  

    }while(np>0);
    inp->end=inp->cur+k; 
    return k;
} 




int next_file(int id)
// moves to next inp file
// uses epoch list
// returns 0 if OK
// 1 if eof
{ 
    double sec;
    int i,e_no;
    struct inpdef *inp=in+cur_inp;
    do{
	e_no=inp->cur_epoch+1;
	if(e_no>inp->start_file[0])
	    return inp->eof=inp_eof=1;
    }while(id>0 && inp->id[e_no]!=id);
    do{
	sec=inp->start_file[e_no];  
	pos_all_files(e_no,sec);
// check if in selection
	e_no=0;
	if(parallel_read)
	    for(i=0;i<no_of_inp;i++)
		e_no+=find_next(in+i,0);
	else
	    e_no=find_next(inp,0);
	if(e_no)
	    return inp_eof=1;
    }while(id>0 && inp->id[e_no]!=id);
    inp->gap=inp->new_epoch=inp->new_out=1;
    return inp_eof=0;
}


int next_selection(int id)
{   
// moves to next selection epoch
// with id (if id==0 all)
// if no selection file 
// calls next file
// returns 0 if OK
// 1 if eof
// uses only first input stream for selections
    struct inpdef *inp=in+cur_inp; 
    
    if(!inp_sel)
	return next_file(id);
    if(read_select_line(inp))
	return inp_eof=inp->eof=1;
//    printf("start=%d\n",inp->e_start_pos);
    inp_eof=inp->eof=-find_next(inp,inp->e_start_pos);
    inp->gap=inp->new_epoch=inp->new_out=1;
    return inp_eof;
}



int next_epoch(int id)
{
// moves to next file or selection  with id 
//dependent on global flag epoch_def

    switch(epoch_def)
    {
	case selections:
	    return next_selection(id);
	case files:
	    return next_file(id);
    }
}
	

void rewind_inp(struct inpdef *inp)
{
    if(inp->sel)rewind(inp->sel);
    inp->e_end_pos=-1;
    inp->end=inp->end_pos=0;
    inp->e_start_pos=-INT_MAX;
    inp->eof=0;
    pos_file(inp,1,inp->start_file[1]);
 //   find_next(inp,0);
}
 
void rewind_all_inp()
{ int i;
 cur_inp=0;
 if(parallel_read)
     for(i=0;i<no_of_inp;i++)
	 rewind_inp(in+i);
 else
     rewind_inp(in);
 inp_eof=0;
}

void new_out_len(struct outdef *out, double win_len, double win_shift)
{
    if(!out || win_shift<=0. )
	print_error(ERROR,"new_out_len","out stream not set or illegal parameters",NULL);
    out->win_len=win_len;
    out->win_shift=win_shift;
    out->win_rate=1.l/win_shift;
}
void set_out_len(struct outdef *out, struct inpdef * inp)
{
// sets out  win and shift same as inp

    if(!out || !inp || !inp->freq || !inp->inp_spoints)
	print_error(ERROR,"set_out_len","inp or out stream not set",NULL);
    out->win_shift=1.l/inp->freq;
    out->win_len=inp->inp_npoints*out->win_shift/inp->inp_spoints;
    out->win_rate=inp->freq;
}
