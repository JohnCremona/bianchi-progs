
# Makefile for bianchi-progs with test programs

# eclib is a requirement.  If installed in the usual place /usr/local
# the following line might not be necessary.  If installed anywhere
# else set the install directory here:

ECLIB_BASE=/usr/local
#ECLIB_BASE=$(HOME)/eclib
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

#  USE_SMATS compiler flag (set/unset in face_relations.h)
#
#  It is recommended set USE_SMATS when finding rational newforms
#  only, since the linear algebra is done correctly modulo p=2^30-35
#  and only lifted to Z when an eigenvector is found. Programs which
#  compute Hecke matrices and their characteristic polynomials such as
#  hecketest.cc try to lift the whole modular symbol space and this
#  easily fails.  In that case unset USE_SMATS, and the linear algebra
#  will be done using multiprecision integer arithmetic instead.  This
#  is much slower: e.g. with level (128), field 1 it takes 20m instead
#  of <1s.

GCC=g++ -std=c++11 -fmax-errors=1
CC = $(GCC)

# to disable checking of assert() use the following:
#OPTFLAG = -DNDEBUG -O3 -fPIC
# to enable checking of assert() use the following:
OPTFLAG = -O3 -Wall -fPIC


# The type of integers used for components of Quad, Qideal, RatQuad
# can be either long or bigint (=NTL's ZZ), these being typedef'd to
# INT, or the class INT (wrapper around FLINT's fmpz_t defined in
# int.h/cc).  Choose which by setting INT_TYPE here.  If the type is
# long then INT_IS_long is defined; this is fastest but results in
# overflow for large levels over larger fields.  If it is ZZ, there is
# no overflow but the code is slower; otherwise the macro FLINT is
# defined, there is no overflow and the code runs about 4 times slower
# than using longs but about 5 times faster than using ZZ.

# Timings for "make check" 22/11/23:
# INT_TYPE=long:  4m 24s =  264s
# INT_TYPE=ZZ:   78m 25s = 4705s (~ 17.8 times slower)
# INT_TYPE=none: 15m 47s =  947s (~  3.6 times slower)

INT_TYPE=
ifeq ($(INT_TYPE), ZZ)
 BASE_TYPE_FLAG = -D INT_IS_ZZ
else
ifeq ($(INT_TYPE), long)
 BASE_TYPE_FLAG = -D INT_IS_long
else
 BASE_TYPE_FLAG = -D FLINT
 FLINT_LDFLAGS = -lflint -lgmp
endif
endif

# NB If used with a multithreaded build of eclib then you MUST define
# USE_BOOST=1 below so that the correct compiler and linker stuff is
# appended below.  Otherwise set USE_BOOST=0 (or do not set it).

#USE_BOOST=1
ifeq ($(USE_BOOST), 1)
 BOOST_ASIO_LIB = -lboost_system-mt
 BOOST_CPPFLAGS =   -DECLIB_MULTITHREAD -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_TR1_UNORDERED_MAP=/\*\*/ -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_UNORDERED_MAP=/\*\*/# -pthread -I/usr/include
 BOOST_SYSTEM_LIB = -lboost_system
 BOOST_THREAD_LIB = -lboost_thread
 BOOST_LDFLAGS = -L/usr/lib -pthread $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
endif

# for profiling:
#CFLAGS = -c -g -pg $(OPTFLAG) $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
#LFLAGS = -pg $(FLINT_LDFLAGS) -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for coverage:
#CFLAGS = -c -g --coverage $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
#LFLAGS = --coverage -fprofile-arcs $(FLINT_LDFLAGS) -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for normal use:
CFLAGS = -c -g $(OPTFLAG) $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
LFLAGS = $(FLINT_LDFLAGS) -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

all: tests

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs0 ccs1 ccs2 ccs3 ccs4 ccs5 ccs6 ccs7
ccs0: int.cc arith_extras.cc intprocs.cc matprocs.cc quads.cc mat22.cc fieldinfo.cc cusp.cc homtest.cc hecketest.cc
ccs1: lf1.cc looper.cc looptest.cc euclid.cc geometry.cc
ccs2: P1N.cc newforms.cc oldforms.cc homspace.cc edge_relations.cc face_relations.cc hecke.cc
ccs3: testlf1.cc makenf.cc pmanin.cc tquads.cc tratquad.cc dimtable.cc dimtabeis.cc dimtabnew.cc dimtabtwist.cc dimtable_all.cc
ccs4: nftest.cc nflist.cc moreap.cc moreap1.cc moreap_loop.cc modularity.cc modularity_modp.cc
ccs5: qideal.cc qidloop.cc primes.cc qidltest.cc qidl_labels.cc
ccs6: hecketest_modp.cc dimtable_modp.cc makenf_modp.cc nflist_modp.cc rewrite_eigs.cc flint_test
ccs7: swan.cc swan_test.cc

headers: arith_extras.h int.h rat.h intprocs.h matprocs.h cusp.h homspace.h lf1.h looper.h P1N.h newforms.h oldforms.h quads.h ratquads.h euclid.h geometry.h qideal.h primes.h qidloop.h mat22.h hecke.h swan.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = fieldinfo tquads qidltest tratquad looptest homtest hecketest makenf moreap moreap1 nftest nflist dimtable dimtable_all dimtabeis dimtabnew dimtabtwist modularity modularity_modp P1Ntest dimtable_modp hecketest_modp makenf_modp makenf_loop nflist_loop rewrite_eigs qidl_labels flint_test swan_test

tests: $(TESTS)

# These are for creation of temporary newforms directories for tests:
DISCS9=3 4 7 8 11 19 43 67 163
DISCSXodd=15 23 31 35 39 47 51 55 59 71 79 83 87 91 95
DISCSXeven=20 24 40 42 52 56 68 84 88 168
DISCS=$(DISCS9) $(DISCSXodd) $(DISCSXeven)

# These control which tests run on which fields:

# All tests: basic arithmetic, homology dimensions, newforms, modularity
FIELDS_full=1 2 3 7 11 19 43 67 163 23 31 47 59 71 79 83 5 14 21 95
#FIELDS_full=

# Basic arithmetic, homology dimensions, newforms
FIELDS_nf=$(FIELDS_full) 6 10 13 15 17 22 35 39 42 51 55 87 91
#FIELDS_nf=$(FIELDS_full)

# Basic arithmetic, homology dimensions
FIELDS_hom=$(FIELDS_nf)
#FIELDS_hom=

# Only basic arithmetic
FIELDS=$(FIELDS_hom)
#FIELDS=

# modtest and symbtest no longer maintained as classes moddata, symbdata are obsolete
BASIC_TESTS = tquads tratquad looptest qidltest
#BASIC_TESTS = tratquad looptest qidltest
HOM_TESTS = homtest dimtable dimtabeis hecketest #dimtable_modp hecketest_modp nflist_modp
NF_TESTS = makenf_loop makenf nftest nflist nflist_loop dimtabnew dimtabtwist moreap moreap1
FULL_TESTS = modularity modularity_modp  #makenf_modp
# global tests are universal, not per field
GLOBAL_TESTS = fieldinfo dimtable_all P1Ntest
ALL_TESTS = $(BASIC_TESTS) $(HOM_TESTS) $(NF_TESTS) $(FULL_TESTS) $(GLOBAL_TESTS)

test_input_dir = testin
test_output_dir = testout

TIMES := $(shell mktemp)

check_run = echo -n "Testing $${prog} for d=$${d}..."; time -o $(TIMES) -f "%Us" ./$${prog} < $(test_input_dir)/$${prog}.$${d}.in > $${prog}.$${d}.out 2>/dev/null && if diff -q $${prog}.$${d}.out $(test_output_dir)/$${prog}.$${d}.out; then echo "$${prog} for d=$${d} completed successfully in " `cat $(TIMES)`;  else echo " ! $${prog} for d=$${d} failed"; diff $${prog}.$${d}.out $(test_output_dir)/$${prog}.$${d}.out; fi || exit $$?


export NF_DIR:=nftmp
check: $(ALL_TESTS)
	 @echo Test setup: create temporary directories
	 rm -f t
	 rm -rf $(NF_DIR)
	 mkdir $(NF_DIR)
	 for d in $(DISCS); do mkdir $(NF_DIR)/2.0.$$d.1; done
	 @echo
	 @echo running global tests...
	 @echo
	 @for d in all; do for prog in $(GLOBAL_TESTS); do $(check_run); done; echo; done
	 @echo
	 @echo running basic tests on fields $(FIELDS)...
	 @echo
	 @for d in $(FIELDS); do for prog in $(BASIC_TESTS); do $(check_run); done; echo; done
	 @echo
	 @echo running basic homspace tests on fields $(FIELDS_hom)...
	 @echo
	 @for d in $(FIELDS_hom); do for prog in $(HOM_TESTS); do $(check_run); done; echo; done
	 @echo
	 @echo running newform tests on fields $(FIELDS_nf)...
	 @echo
	 @for d in $(FIELDS_nf); do for prog in $(NF_TESTS); do $(check_run); done; echo; done
	 @echo
	 @echo running full tests on fields $(FIELDS_full)...
	 @echo
	 @for d in $(FIELDS_full); do for prog in $(FULL_TESTS); do $(check_run); done; echo; done
	 @echo
	 @echo Tidy up: remove temporary directories and output test files
	 rm -rf $(NF_DIR)
	 rm -f *.out
	 @echo Tests completed

clean:
	rm -f $(TESTS)
	rm -f *.o *~ *.testout

OBJS = int.o arith_extras.o intprocs.o quads.o matprocs.o euclid.o geometry.o looper.o homspace.o \
       newforms.o oldforms.o edge_relations.o face_relations.o hecke.o qideal.o qidloop.o \
       primes.o mat22.o ratquads.o cusp.o P1N.o swan.o

objs: $(OBJS)

makenf_loop.o: makenf.cc
	$(CC) -DLOOPER $(CFLAGS) makenf.cc -o makenf_loop.o

nflist_loop.o: nflist.cc
	$(CC) -DLOOPER $(CFLAGS) nflist.cc -o nflist_loop.o

flint_test.o: flint_test.cc int.cc int.h rat.h
	$(CC) $(CFLAGS) $<

swan_test.o: swan_test.cc int.cc int.h rat.h swan.o
	$(CC) $(CFLAGS) $<

tquads: tquads.o $(OBJS)
	$(CC) -o tquads tquads.o $(OBJS) $(LFLAGS)

P1Ntest: P1Ntest.o $(OBJS)
	$(CC) -o P1Ntest P1Ntest.o $(OBJS) $(LFLAGS)

fieldinfo: fieldinfo.o $(OBJS)
	$(CC) -o fieldinfo fieldinfo.o $(OBJS) $(LFLAGS)

makenf: makenf.o $(OBJS)
	$(CC) -g -o makenf makenf.o $(OBJS) $(LFLAGS)

makenf_modp: makenf_modp.o $(OBJS)
	$(CC) -g -o makenf_modp makenf_modp.o $(OBJS) $(LFLAGS)

makenf_loop: makenf_loop.o $(OBJS)
	$(CC) -g -o makenf_loop makenf_loop.o $(OBJS) $(LFLAGS)

pmanin: pmanin.o $(OBJS)
	$(CC) -o pmanin pmanin.o $(OBJS) $(LFLAGS)

testlf1: testlf1.o lf1.o $(OBJS)
	$(CC) -o testlf1 testlf1.o lf1.o $(OBJS) $(LFLAGS)

nftest: nftest.o $(OBJS)
	$(CC) -o nftest nftest.o $(OBJS) $(LFLAGS)

nflist: nflist.o $(OBJS)
	$(CC) -o nflist nflist.o $(OBJS) $(LFLAGS)

nflist_loop: nflist_loop.o $(OBJS)
	$(CC) -o nflist_loop nflist_loop.o $(OBJS) $(LFLAGS)

nflist_modp: nflist_modp.o $(OBJS)
	$(CC) -o nflist_modp nflist_modp.o $(OBJS) $(LFLAGS)

moreap: moreap.o $(OBJS)
	$(CC) -o moreap moreap.o $(OBJS) $(LFLAGS)

moreap1: moreap1.o $(OBJS)
	$(CC) -o moreap1 moreap1.o $(OBJS) $(LFLAGS)

moreap_loop: moreap_loop.o $(OBJS)
	$(CC) -o moreap_loop moreap_loop.o $(OBJS) $(LFLAGS)

modularity: modularity.o $(OBJS)
	$(CC) -o modularity modularity.o $(OBJS) $(LFLAGS)

modularity_modp: modularity_modp.o $(OBJS)
	$(CC) -o modularity_modp modularity_modp.o $(OBJS) $(LFLAGS)

looptest: looptest.o $(OBJS)
	$(CC) -o looptest looptest.o $(OBJS) $(LFLAGS)

tratquad: tratquad.o $(OBJS)
	$(CC) -o tratquad tratquad.o $(OBJS) $(LFLAGS)

homtest: homtest.o $(OBJS)
	$(CC) -o homtest homtest.o $(OBJS) $(LFLAGS)

dimtable_modp: dimtable_modp.o $(OBJS)
	$(CC) -o dimtable_modp dimtable_modp.o $(OBJS) $(LFLAGS)

dimtable: dimtable.o $(OBJS)
	$(CC) -o dimtable dimtable.o $(OBJS) $(LFLAGS)

dimtable_all: dimtable_all.o $(OBJS)
	$(CC) -o dimtable_all dimtable_all.o $(OBJS) $(LFLAGS)

dimtabeis: dimtabeis.o $(OBJS)
	$(CC) -o dimtabeis dimtabeis.o $(OBJS) $(LFLAGS)

dimtabnew: dimtabnew.o $(OBJS)
	$(CC) -o dimtabnew dimtabnew.o $(OBJS) $(LFLAGS)

dimtabtwist: dimtabtwist.o $(OBJS)
	$(CC) -o dimtabtwist dimtabtwist.o $(OBJS) $(LFLAGS)

hecketest: hecketest.o $(OBJS)
	$(CC) -o hecketest hecketest.o $(OBJS) $(LFLAGS)

hecketest_modp: hecketest_modp.o $(OBJS)
	$(CC) -o hecketest_modp hecketest_modp.o $(OBJS) $(LFLAGS)

roundtest: roundtest.o quads.o
	$(CC) -o roundtest roundtest.o quads.o $(LFLAGS)

qidltest: qidltest.o $(OBJS)
	$(CC) -o qidltest qidltest.o $(OBJS) $(LFLAGS)

qidl_labels: qidl_labels.o $(OBJS)
	$(CC) -o qidl_labels qidl_labels.o $(OBJS) $(LFLAGS)

rewrite_eigs: rewrite_eigs.o $(OBJS)
	$(CC) -o rewrite_eigs rewrite_eigs.o $(OBJS) $(LFLAGS)

flint_test: flint_test.o int.h rat.h
	$(CC) -o flint_test flint_test.o int.o $(LFLAGS)

swan_test: swan_test.o $(OBJS) rat.h
	$(CC) -o swan_test swan_test.o $(OBJS) $(LFLAGS)

# DEPENDENCIES
#
# recreate with
# for f in *.cc; do g++ -MM -std=c++11 ${f}; done > Makefile.deps
#
include Makefile.deps

# Some tables

# 1. Data tabulated in 1984 Compositio paper Tables 3.{2,3,4,5,6}.1:

paperdims: paperdims1 paperdims2 paperdims3 paperdims7 paperdims11
paperdims1: dimtable
	echo 1 0 0 1 500 | ./dimtable | awk '$$5>0' | tail -n +3 > paperdims.1.out
paperdims2: dimtable
	echo 2 0 0 1 300 | ./dimtable | awk '$$5>0' | tail -n +3 > paperdims.2.out
paperdims3: dimtable
	echo 3 0 0 1 500 | ./dimtable | awk '$$5>0' | tail -n +3 > paperdims.3.out
paperdims7: dimtable
	echo 7 0 0 1 200 | ./dimtable | awk '$$5>0' | tail -n +3 > paperdims.7.out
paperdims11: dimtable
	echo 11 0 0 1 200 | ./dimtable | awk '$$5>0' | tail -n +3 > paperdims.11.out
