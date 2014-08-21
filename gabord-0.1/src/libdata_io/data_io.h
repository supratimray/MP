/*************************************************
*    Copyright Piotr J. Franaszczuk 2005         *
*    ERL JHU
**************************************************/
#ifndef _DATA_IO_H_

#define _DATA_IO_H_
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#ifndef WINDOWS
	#include <unistd.h>
#endif
#include <fcntl.h>
#include "cfg_io.h"

// MEM macros require sname defined
#define DMEM(var,n)  if(!(var=(double*)calloc((n),sizeof(double)))) mem_error(sname,"var");
#define IMEM(var,n)  if(!(var=(int*)calloc((n),sizeof(int)))) mem_error(sname,"var");
#define SMEM(var,n)  if(!(var=(short*)calloc((n),sizeof(short)))) mem_error(sname,"var");
#define MEM(type,var,n)  if(!(var=(type*)calloc((n),sizeof(type)))) mem_error(sname,"var");
//#define NAME_LEN 50

#define NO_MTG  0  // no montage change all channels included
#define MTG_SEL 1  // no reference change but only some channels selected
#define MTG_REF 2  // reference change in some channels
#define MTG_AVR 3  // some channels have average reference

enum e_fmt {little_end,big_end,
	    text,
	    bin_char,bin_byte,                   // 8 bytes
	    bin_int12,                           // 12 bytes packed            
	    bin_short,bin_ushort,bin_edf,       //  16 bytes
	    bin_int,bin_uint,bin_float,         // 32 bytes
	    bin_double,                         // 64 bytes
	    mode_serial,mode_parallel,
	    selections,files
};


#if(__BYTE_ORDER == 1234 || _WIN32)
#define MY_BYTE_ORDER little_end   //Intel
#else
#define MY_BYTE_ORDER big_end      //MIPS
#endif


#ifndef uint
#define uint unsigned int
#endif

#ifndef ushort
#define ushort unsigned short
#endif

struct mont {
  int type;  // type of montage
  int nchans; // number of resulting channels for analysis
  int mchans; // number of channels needed for computing resulting channels
  char **chan_labels; // pointer to array of nchans labels [nchans]
  char **chan_units;  // pointer to array of nchans labels for units
  //  char **var_labels; // pointer to array of variable labels [nvars]
  //  char **var_units;  // pointer to array of variable labels for units
  double *chan_calib; //pointer to array of calibration factors (if present) 
  int *chan_no; // array of indexes to channels in array of selected (dbuf) [nchans]
  int *ref_no;  // references indexes to  needed channels[nchans]
  int *ch_no;  // indexes to original channels (rbuf)  [mchans]
  int navg;   // max number of averages used i.e. largest N in AVGN
  int **lst_avg; // arrays of channels used for average [navg][avg_no] (dbuf)
  int *avg_no; // no of channels averaged [navg]
  char **avg_defs;  // labels of channels used in def of AVG
  // lists of channels in avg [navg][*] if avg_lst[I]>0 compute AVGI
  // AVG0 and AVG1 predefined, but if not used not computed
    
};


struct hfmt {
  char *title;            
  double samp_rate;
  int   bytes_per_win;
  int   nchans;
  int   byte_order;
  int   file_fmt;
  int   offset;
  int   chans_per_file;    // no_of chans per file
  double win_len;
  double win_shift;
  double win_rate;         //=1./ win_shift (it helps with rounding errors )    
  int   no_of_vars;        // no of vars/chan/win  0 if variable
  int   bytes_per_var;     // determined by file_fmt unless text
  int   fil_ch_no;         // start no of ext to name if chans_per_file==1
  char* lst_name;          // name of list file 
  char* templ;             // name  template
  char ** chan_labels;
  char **chan_units; 
  char **avg_defs;
  double * chan_calib;
  char* mtg_name, *vmtg_name;
 
};
  
struct pinfo { char * ID, *fname, *lname, *dob ;};

struct inpdef {
  char *title;
  FILE *fil,*lst,*sel,**files;
  char *path,*hdr_name,*lst_name,*sel_name,*mtg_name,*vmtg_name;
  int pipe;            // =1 if fifo pipe
  int new_epoch;       // flag indicating new epoch started in last find
  int var_len;         // variable length records
  int (*skip_fun)(struct inpdef *inp, int pos);
  int (*read_fun)(struct inpdef *inp, char*sbuf,int n);
  int gap;            // ==1 if gap in files or new epoch in last find
  int new_out;        // ==1 if write should start new out file
  int   file_fmt;
  void (*swap_fun)(short* b,int n,int m); 
  // function to swap if byte order not native
  //  (also decode_tf)
  void (*select_fun)(struct inpdef *inp,char*sbuf,double *db,int n); 
  // to put selected channels in dbuf
  void *rbuf;          // address of input buffer
  double *dbuf;        // pointer to double buffer with selected channels (and averages)
  short *edf_buf,*edf_sbuf;
  int edf_len,rec_len;    // for edf input files 
  double f_start,f_end,f_cur,e_start,e_end; // start and end of current file or epoch in secs
  int pos,cur,end,end_pos,e_end_pos,e_start_pos;  // positions relative to f_start in inp windows
  int eof;             // eof stream == 1 when it is end of list or selections
  char  name[NAME_LEN];
  int nchans;          // nchans in INPUT file if multiple files this is only needed chans i.e. =mtg.mchans
  int no_of_chans;     // no of chans in DATA window same as mtg.nchans
  double samp_rate;     // original samp rate
  int inp_npoints;     //  it is always == INPUT npoints(hmft.win_len*freq)
  int inp_spoints;     //  it is always == INPUT spoints(hmft.spoints)
  double freq;         // effective freq for INPUT = 1/(INPUT win_shift) 
  int   bytes_per_win; // number of bytes per INPUT window do not use after read
  int   bytes_per_var; // determined by file_fmt unless text
  int   offset;        // offset of data in file (bytes)cc -o sumemp -O3 sumemp.c data_io.c cfg_io.c -lgsl -lgslcblas -lblas -lg2c -lm
  int zeros;            // if !=0 put zeros for names starting with # in list file, otherwise skip those files(default) if >0 put zeros in current
  char * epoch;        // epoch for analysis in same format as select file )ignored if select file present) 
  int   no_of_vars;    // no of vars/chan/win  
                       // no of vars/per INPUT/DATA window = 
                       // bytes_per_win/nchans/sizeof(fmt) 
  int*   no_chan_vars;    // table for no_of_vars for each channel 
  int   win_size;      // size of window for analysis in vars 
                       // (=no_of_chans*no_of_vars*npoints)
  int   npoints;       // same as win_len but in points
  int   spoints;       // same as win_shift but in points
  int rate_flag;       //   win_rate takes precedence over win_shift     
  double win_len;       // len of DATA window in sec
  double win_shift;     // shift for next DATA window in sec
  double win_rate;      // =1/win_shift
  int   bytes_per_chan; 
  int   bytes_win_size;// size of window after unpacking 
  int   dwin_size;     // =no_of_vars*(mtg.mchans+1+mtg.navg) // if only selection mtg (there is no +1)
  int   calib;         // if 0 no calibration even if exist in header
  int   fil_ch_no;     // start no of ext to name if chans_per_file==-1
  int   *ch_no;        // original channel numbers copy of mtg.ch_no if chans_per_file==-1, otherwise NULL
  int  nchvars;       // no_of_chans*no_of_vars
  char **labels;     // data labels including chans and vars
  char **units;     //data labels including chans and vars  

  double *list,*start_file; // lists for direct access reading (no sequential)
  int cur_epoch;  // for direct access counted from 0
  char ** fname;
    int id_select;     // id to select from select or list line
  int sel_id;          // id  from selection line
  int sel_out_id;     // the previous value needed for output
  char *sel_out_text;
  char* sel_text;        // string from selection line
  int * id;            // array of id from list file for direct access reading
  char **text;         // array for strings from list lines 
  int sel_comment;     // if 1 sel_comment is taken from select file,  otherwise from list file

  struct mont mtg;     // structure with all information about chans in montage
  struct mont vmtg;    // structure with all information about vars in montage
  struct pinfo p_info;
 
};
// bytes_per_win=bytes_per_chan*nchans=bytes_per_var*no_of_vars*nchans

struct outdef { 
  char *title;
  FILE *fil,*lst,*hdr,**files;
  int pipe;
  int var_len;           // variable no  of variables here max
  char name[NAME_LEN];
  char *path,*hdr_name,*lst_name;
  double start,len;
  char *templ,*ext;
  FILE* alarm;           // this is file for writing alarm file (dflt for int12)
  int   bytes_per_win;   // number of bytes per window =nchans*npoints*b_p_var 
  int   no_of_vars;      // number of vars per chan if  variable output as negative
  int*   no_chan_vars;    // table for no_of_vars for each channel 
  int   bytes_per_var;   // determined by file_fmt unless text
  uint   bytes;           // bytes counter
  uint*  mbytes;
  int   byte_order;
  int   file_fmt;        // if not text it may also determine bytes_per_unit 
  int   no_of_chans;     // no of output channels (dflt in->no_of_chans )
  int   start_chan;      // start channel for output (dflt 1 (first))
  int   all_chans;     // no of chans in all out_streams to split  (dflt global all_chans)
  int   chans_per_file;  // no_of chans per file
  int   file_no;         // current file_no;
  double samp_rate;       // sampling rate for output header usually same as input
  uint   max_len;         // max len of output file in output windows
  int   new_out;         // 1 if write opened new file
  int   count;           // counter for no of variables written
  int join;              // if not zero makes continous file even if there are gaps, use only for export!!!
  void (*swap_fun)(short* b,int n,int m);
  // function to swap if byte order not native
  void (*write_data)(struct outdef *out,double *buf,int npoints); // writes data in correct format
  void (*write_chan_data)(struct outdef *out,double *buf,int nv,int chan,int nvars ); 
  double win_len;         // len of OUTPUT window in sec
  double win_shift;       // shift for next OUTPUT window in sec
  double  win_rate;       // win_rate =1./win_shift
  double offset;          // offset in start of out file (in sec.)  dfl =0. Use for filter etc...
  int   rate_flag;        // if ==1 win_rate takes precedence over win_shift
  int   inp_no;          // index to inp stream associated with this output
  int   calib;           // if 1 copy calibration from input mtg dflt 0
  char **chan_units; 
  char **chan_labels;
  char **var_units; 
  char **var_labels;
  int fil_ch_no;         // start no of ext to name if chans_per_file==1
  int templ_dig;           // no of digits in file number in template (dflt 4)
  int sel_id;          // id to database record from selction line
  char* sel_text;        // string from selection line
};



#ifdef CPP
#define _DECL extern "C"
#else
#define _DECL extern
#endif

_DECL int prec(double x);
_DECL  struct inpdef * get_inp_str(struct inpdef* inp,char* title);
_DECL  struct outdef * get_out_str(struct outdef *out,char* title);
_DECL void concat(struct inpdef *inp);
_DECL void write_out_hdr(struct outdef *out);
_DECL int read_mtg(char *name, struct mont*, struct hfmt *); 
_DECL int read_cfg(char *name,struct section *head);

_DECL void swap2(short *b, int n, int m);
_DECL void swap4(short *b, int n, int m);
_DECL void swap8(short *b, int n, int m);
_DECL void decode_tf(short *b, int n, int m);
_DECL void decode_tf_swap(short *b, int n, int m);



_DECL void compute_avg(struct inpdef *inp,double *db,int nn);
_DECL void select_chans_short_calib(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_ushort_calib(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_short(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_ushort(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_byte(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_char(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_int(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_uint(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_float(struct inpdef *inp,char*sbuf,double* db,int n);
_DECL void select_chans_double(struct inpdef *inp,char*sbuf,double* db,int n);

_DECL double* get_epochs(struct inpdef *inp,int merge_sel);
_DECL int find_in_list(struct inpdef *inp,double sec,double len);
_DECL int find_in_list_all(double,double);
_DECL int pos_file(struct inpdef *inp,int j, double sec);
_DECL int pos_all_files(int,double);
_DECL int find(struct inpdef *inp,double sec);
_DECL int find_next(struct inpdef *inp,int n);
_DECL int next_epoch(int);
_DECL int next_selection(int);
_DECL int next_file(int);
_DECL int find_forward(struct inpdef *inp,int n);
_DECL double get_time(char *time, double tod);
_DECL int    read_list_line(struct inpdef *inp);
_DECL int    read_select_line(struct inpdef *inp);
_DECL int    write_select_line(FILE *sel,double dt,double t_cur,double len,int no, char* type);
_DECL int    write_list_line(struct outdef *out);
_DECL void   write_out_hdr(struct outdef *out);
_DECL void   new_out_file(struct outdef *out);
_DECL void   new_out_file_chan(struct outdef *out,int chan);
_DECL int    new_inp_file(struct inpdef *inp,int flag);

_DECL int    init_io(int argc,char **argv, struct section *list, int flag);
_DECL void   close_io();
_DECL void   write_chan_data(struct outdef *out,void *buf,int nv,int chan,int nvars);
_DECL void   write_data_raw(struct outdef *out,void *buf,int npoints,int nvars,int size);
_DECL int    read_next_win(struct inpdef *inp,int n,int part);


_DECL int read_win_fix(struct  inpdef *inp, char* sbuf, int n);
_DECL int read_win_edf(struct  inpdef *inp, char* sbuf, int n);
_DECL int read_win_fix_chan(struct  inpdef *inp, char* sbuf, int n,int chan);
_DECL int read_win_fix_many(struct inpdef *inp,char* sbuf,int n);
_DECL int read_win_var(struct  inpdef *inp, char* sbuf, int n);
_DECL int read_win_var_chan(struct  inpdef *inp, char* sbuf, int n,int chan);
_DECL int read_win_var_many(struct inpdef *inp,char* sbuf,int n);

_DECL int skip_fix(struct inpdef *inp,int pos);
_DECL int skip_edf(struct inpdef *inp,int pos);
_DECL int skip_fix_many(struct inpdef *inp,int  pos);
_DECL int skip_var(struct inpdef *inp,int  pos);
_DECL int skip_var_many(struct inpdef *inp,int pos);
_DECL void get_vars(struct inpdef *inp, void *buf,int n, int ch, int var, int nvars);
_DECL int get_double_chan(struct inpdef *inp,double *data,int n,int chan);
_DECL int get_float_chan(struct inpdef *inp,float *data,int n,int chan);
_DECL int get_int_chan(struct inpdef *inp,int *data,int n,int chan);
_DECL int get_short_chan(struct inpdef *inp,short *data,int n,int chan);
_DECL int get_uint_chan(struct inpdef *inp,uint *data,int n,int chan);
_DECL int get_ushort_chan(struct inpdef *inp,ushort *data,int n,int chan);
_DECL int get_char_chan(struct inpdef *inp,char *data,int n,int chan);
_DECL int get_byte_chan(struct inpdef *inp,unsigned char *data,int n,int chan);
_DECL int get_raw_chan(struct inpdef *inp,void *data,int n,int chan);

_DECL int get_double_data_mtg(struct inpdef *inp, double *data,int n,struct mont *mtg);
_DECL int get_double_data(struct inpdef *inp, double *data,int n);
_DECL int get_float_data(struct inpdef *inp, float *data,int n);
_DECL int get_int_data(struct inpdef *inp, int *data,int n);
_DECL int get_short_data(struct inpdef *inp, short *data,int n);
_DECL int get_uint_data(struct inpdef *inp,uint *data,int n);
_DECL int get_ushort_data(struct inpdef *inp,ushort *data,int n);
_DECL int get_char_data(struct inpdef *inp,char *data,int n);
_DECL int get_byte_data(struct inpdef *inp,unsigned char *data,int n);
_DECL int get_raw_data(struct inpdef *inp,void *data,int n);

_DECL int get_double_all_data( double *data,int n);
_DECL int get_float_all_data( float *data,int n);
_DECL int get_int_all_data( int *data,int n);
_DECL int get_short_all_data( short *data,int n);
_DECL int get_uint_all_data(uint *data,int n);
_DECL int get_ushort_all_data(ushort *data,int n);
_DECL int get_char_all_data(char *data,int n);
_DECL int get_byte_all_data(unsigned char *data,int n);
_DECL int get_raw_all_data(void *data,int n);

_DECL int read_double_data(struct inpdef *inp, double *data,int n,int part);
_DECL int read_float_data(struct inpdef *inp, float *data,int n,int part);
_DECL int read_int_data(struct inpdef *inp, int *data,int n,int part);
_DECL int read_short_data(struct inpdef *inp, short *data,int n,int part);
_DECL int read_uint_data(struct inpdef *inp,uint *data,int n,int part);
_DECL int read_ushort_data(struct inpdef *inp,ushort *data,int n,int part);
_DECL int read_char_data(struct inpdef *inp,char *data,int n,int part);
_DECL int read_byte_data(struct inpdef *inp,unsigned char *data,int n,int part);

_DECL int read_double_all_data(double *data,int n,int part);
_DECL int read_float_all_data(float *data,int n,int part);
_DECL int read_int_all_data(int *data,int n,int part);
_DECL int read_short_all_data(short *data,int n,int part);
_DECL int read_uint_all_data(uint *data,int n,int part);
_DECL int read_ushort_all_data(ushort *data,int n,int part);
_DECL int read_char_all_data(char *data,int n,int part);
_DECL int read_byte_all_data(unsigned char *data,int n,int part);
_DECL int read_raw_all_data(void *data,int n,int part);

_DECL void   write_chan_double_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_float_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_int_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_uint_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_short_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_ushort_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_char_data(struct outdef *out,double *buf,int len,int chan,int nvars);
_DECL void   write_chan_byte_data(struct outdef *out,double *buf,int len,int chan,int nvars);

_DECL void   write_double_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_float_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_int_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_uint_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_short_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_ushort_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_char_data(struct outdef *out,double *buf,int npoints);
_DECL void   write_byte_data(struct outdef *out,double *buf,int npoints);
#define      set_new_out(out) in[out->inp_no].new_out=1;

_DECL void new_inp_size(struct inpdef *in, int n);
_DECL void set_out_fmt(struct outdef *out, int nchans, int nvars, int fmt, int files);
_DECL void write_all_data(double *buf,int npoints);
_DECL void write_all_data_raw(void *buf,int npoints,int nvars,int varsize);

_DECL int get_selection(FILE* fil,double f,double *beg, double *end, int*type_no,char **type);

_DECL void new_out_len(struct outdef*,double,double);
_DECL void set_out_len(struct outdef*,struct inpdef*);
#ifdef _WIN32
#define rint(x) floor(x+0.5)
#define bzero(b,n) memset(b,0,n)
#endif
#endif
