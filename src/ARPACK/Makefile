############################################################################
#
#  Program:         ARPACK
#
#  Module:          Makefile
#
#  Purpose:         Sources Makefile
#
#  Creation date:   February 22, 1996
#
#  Modified:        September 9, 1996
#
#  Send bug reports, comments or suggestions to arpack.caam.rice.edu
#
############################################################################
#\SCCS Information: @(#)
# FILE: Makefile   SID: 2.1   DATE OF SID: 9/9/96   RELEASE: 2

#include ./ARmake.inc
RANLIB = $(AR)
ARPACKLIB = libarpack.a
ARFLAGS = rcs

############################################################################
#  To create or add to the library, enter make followed by one or
#  more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is defined by $(ARPACKLIB) in
#  ../ARmake.inc and is created at the next higher directory level.

OBJS = icnteq.o icopy.o iset.o iswap.o ivout.o second_NONE.o

SOBJ  = sgetv0.o slaqrb.o sstqrb.o ssortc.o ssortr.o sstatn.o sstats.o \
	snaitr.o snapps.o snaup2.o snaupd.o snconv.o sneigh.o sngets.o \
	ssaitr.o ssapps.o ssaup2.o ssaupd.o ssconv.o sseigt.o ssgets.o \
	sneupd.o sseupd.o ssesrt.o svout.o smout.o

DOBJ  = dgetv0.o dlaqrb.o dstqrb.o dsortc.o dsortr.o dstatn.o dstats.o \
	dnaitr.o dnapps.o dnaup2.o dnaupd.o dnconv.o dneigh.o dngets.o \
	dsaitr.o dsapps.o dsaup2.o dsaupd.o dsconv.o dseigt.o dsgets.o \
	dneupd.o dseupd.o dsesrt.o dvout.o dmout.o

COBJ  = cnaitr.o cnapps.o cnaup2.o cnaupd.o cneigh.o cneupd.o cngets.o \
        cgetv0.o csortc.o cstatn.o cvout.o cmout.o

ZOBJ  = znaitr.o znapps.o znaup2.o znaupd.o zneigh.o zneupd.o zngets.o \
        zgetv0.o zsortc.o zstatn.o zvout.o zmout.o

.f.o:
	$(F77) $(FFLAGS) -c $<

all: single double complex complex16

single: $(SOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(SOBJ) $(OBJS)

double: $(DOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(DOBJ) $(OBJS)

complex: $(COBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(COBJ) $(OBJS)

complex16: $(ZOBJ) $(OBJS)
	$(AR) $(ARFLAGS) $(ARPACKLIB) $(ZOBJ) $(OBJS)
#
#  clean	- remove all object files
#
clean:
	$(RM) -f *.o
	$(RM) -f $(ARPACKLIB)

