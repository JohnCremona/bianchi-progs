#Makefile for QUADS directory with test programs

METHOD = 2

GCC=g++
CC = $(GCC)
OPTFLAG = -g -O2 -Wall -DMETHOD=$(METHOD)

INCDIR = ../include
LIBDIR = ../lib

IOLIB = -lstdc++  # for gcc-2.7.2 and above

LIDIALIBDIR = /usr/local/LiDIA/lib
LIDIAINCDIR = /usr/local/LiDIA/include

NTLINCDIR = $(SAGE_ROOT)/local/include
NTLLIBDIR = $(SAGE_ROOT)/local/lib

PARIINCDIR = $(SAGE_ROOT)/local/include
PARILIBDIR = $(SAGE_ROOT)/local/lib

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

#
# possible values for ARITH are
#  
#  LIDIA
#  LIDIA_INTS
#  NTL
#  NTL_INTS
#
# override default with "make ARITH=LIDIA" or "make ARITH=LIDIA_INTS"
# or set ARITH as a shell environment variable
#
#
ifndef ARITH
	ARITH=NTL_INTS
endif

%:: RCS/%,v
#  commands to checkout RCS files.  This is the default:
#       test -f $@ || $(CO) $(COFLAGS) $< $@
# and here is my version
	$(CO) $(COFLAGS) $< $@


include targets

shar: sources
	rm -f sharfile
	shar -i files -o sharfile
	compress sharfile

