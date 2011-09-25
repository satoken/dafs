SUBDIRS = contrafold probconsRNA contralign partalign
bin_PROGRAMS = dusaf

dusaf_SOURCES = dusaf.cpp typedefs.h fa.cpp fa.h align.h align.cpp fold.h fold.cpp new_param.c new_param.h ip.cpp ip.h
dusaf_LDADD = contrafold/libcontrafold.a contralign/libcontralign.a probconsRNA/libprobcons.a partalign/libpartalign.a
