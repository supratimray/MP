bin_PROGRAMS = gabord


# set the include path found by configure
INCLUDES = -I$(top_srcdir)/src/libdata_io $(all_includes)

# the library search path.
#gabord_LDFLAGS = $(all_libraries) -R /usr/local/lib
gabord_SOURCES = addons.c complex_op.c cwt1d.h fft.c gabord.c gb_bkold.c \
	gb_buildbook.c gb_corr.c gb_decomp.c gb_filter.c gb_oper.c math.c mpp.h sig_alloc.c \
	sig_fct10.c signals.h struct_alloc.c structpath.h update.c
#gabord_LDADD = $(top_builddir)/src/libdata_io/libdata_io.la -lm
SUBDIRS = libdata_io
gabord_LDADD = $(top_builddir)/src/libdata_io/libdata_io.a -lm
