SUBDIRS = contrafold probconsRNA contralign
bin_PROGRAMS = dusaf

dusaf_SOURCES = dusaf.cpp typedefs.h fa.cpp fa.h \
	align.h align.cpp needleman_wunsch.h needleman_wunsch.cpp \
	fold.h fold.cpp nussinov.h nussinov.cpp ipknot.h ipknot.cpp \
	alifold.cpp alifold.h boltzmann_param.c boltzmann_param.h \
	ip.cpp ip.h cmdline.c cmdline.h
dusaf_LDADD = contrafold/libcontrafold.a contralign/libcontralign.a probconsRNA/libprobcons.a
EXTRA_DIST = gnugetopt.h getopt.c getopt1.c

dusaf.cpp: cmdline.h
cmdline.c cmdline.h: dusaf.ggo
	gengetopt --unamed-opts < dusaf.ggo
