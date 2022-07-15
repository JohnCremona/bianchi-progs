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
#OPTFLAG = -DNDEBUG -O3 -Wall -fPIC
# to enable checking of assert() use the following:
OPTFLAG = -O3 -Wall -fPIC


# The type of integers used for components of Quad, Qideal, RatQuad
# can be either long or bigint (=NTL's ZZ), and is typedef'd to QUINT
# in the code.  By default the type is long and QUINT_IS_long is
# defined; this can result in overflow for large levels over larger
# fields.  Change this to 1 (and make clean and rebuild) to compile
# using ZZ as base integer type for Quads instead.  That results in
# slower code, but it does not overflow!
BASE_TYPE_ZZ=0
ifeq ($(BASE_TYPE_ZZ), 1)
 BASE_TYPE_FLAG = -D QUINT_IS_ZZ
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
#LFLAGS = -pg -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for normal use:
CFLAGS = -c -g $(OPTFLAG) $(BOOST_CPPFLAGS) $(BASE_TYPE_FLAG) -I$(INCDIR)
LFLAGS = -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

all: tests

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs0 ccs1 ccs2 ccs3 ccs4 ccs5
ccs0: intprocs.cc matprocs.cc quads.cc mat22.cc fieldinfo.cc cusp.cc homtest.cc hecketest.cc
ccs1: lf1.cc looper.cc looptest.cc euclid.cc geometry.cc
ccs2: P1N.cc newforms.cc oldforms.cc homspace.cc edge_relations.cc face_relations.cc hecke.cc
ccs3: symb.cc testlf1.cc makenf.cc pmanin.cc tmquads.cc tquads.cc tratquad.cc dimtable.cc dimtabeis.cc dimtabnew.cc dimtabtwist.cc
ccs4: nftest.cc nflist.cc moreap.cc moreap1.cc moreap_loop.cc modularity.cc modularity_modp.cc
ccs5: qideal.cc qidloop.cc primes.cc qidltest.cc hecketest_modp.cc dimtable_modp.cc makenf_modp.cc nflist_modp

headers: intprocs.h matprocs.h cusp.h homspace.h lf1.h looper.h P1N.h mquads.h newforms.h oldforms.h quads.h ratquads.h symb.h euclid.h geometry.h qideal.h primes.h qidloop.h mat22.h hecke.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = fieldinfo tquads qidltest tratquad looptest homtest hecketest makenf moreap moreap1 nftest nflist dimtable dimtabeis dimtabnew dimtabtwist modularity modularity_modp P1Ntest dimtable_modp hecketest_modp makenf_modp makenf_loop nflist_loop
tests: $(TESTS)

# These are for creation of temporary newforms directories for tests:
DISCS9=3 4 7 8 11 19 43 67 163
DISCSXodd=15 23 31 35 39 47 51 55 59 71 79 83 87 91 95
DISCSXeven=20 24 40 42 52 56 68 84 88
DISCS=$(DISCS9) $(DISCSXodd) $(DISCSXeven)

# These control which tests run on which fields:

# All tests: basic arithmetic, homology dimensions, newforms, modularity
FIELDS_full=1 2 3 7 11 19 43 67 163 23 31 47 59 71 79 83
FIELDS_full=

# Basic arithmetic, homology dimensions, newforms
FIELDS_nf=5 6 10 13 14 15 17 22
#FIELDS_nf=5 6 10 13 14 15 17

# Basic arithmetic, homology dimensions
FIELDS_hom=21 35 39 42 51 55 87 91 95
#FIELDS_hom=

# Only basic arithmetic
FIELDSX=
#FIELDSX=

FIELDS=$(FIELDS_full) $(FIELDS_hom) $(FIELDSX)
#FIELDS=

# modtest and symbtest no longer maintained as classes moddata, symbdata are obsolete
BASIC_TESTS = fieldinfo tquads tratquad looptest P1Ntest qidltest
HOM_TESTS = homtest dimtable dimtabeis hecketest #dimtable_modp hecketest_modp nflist_modp
NF_TESTS = $(HOM_TESTS) makenf makenf_loop nftest nflist nflist_loop dimtabnew
FULL_TESTS = $(NF_TESTS) moreap moreap1 modularity modularity_modp  #makenf_modp
ALL_TESTS = $(BASIC_TESTS) $(FULL_TESTS)

test_input_dir = testin
test_output_dir = testout

check_run = echo -n "Testing $${prog} for d=$${d}..."; time -o /dev/tty -f "runtime was %Us..." ./$${prog} < $(test_input_dir)/$${prog}.$${d}.in > $${prog}.$${d}.out 2>/dev/null && if diff -q $${prog}.$${d}.out $(test_output_dir)/$${prog}.$${d}.out; then echo " - $${prog} for d=$${d} completed successfully";  else echo " ! $${prog} for d=$${d} failed"; diff $${prog}.$${d}.out $(test_output_dir)/$${prog}.$${d}.out; fi || exit $$?


export NF_DIR:=nftmp
check: $(ALL_TESTS)
	 @echo Test setup: create temporary directories
	 rm -f t
	 rm -rf $(NF_DIR)
	 mkdir $(NF_DIR)
	 for d in $(DISCS); do mkdir $(NF_DIR)/2.0.$$d.1; done
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

OBJS = quads.o intprocs.o matprocs.o euclid.o geometry.o looper.o homspace.o \
       newforms.o oldforms.o edge_relations.o face_relations.o hecke.o qideal.o qidloop.o \
       primes.o mat22.o ratquads.o cusp.o P1N.o

tquads: tquads.o $(OBJS)
	$(CC) -o tquads tquads.o $(OBJS) $(LFLAGS)

P1Ntest: P1Ntest.o $(OBJS)
	$(CC) -o P1Ntest P1Ntest.o $(OBJS) $(LFLAGS)

fieldinfo: fieldinfo.o $(OBJS)
	$(CC) -o fieldinfo fieldinfo.o $(OBJS) $(LFLAGS)

tmquads: tmquads.o mquads.o
	$(CC)  -o tmquads tmquads.o mquads.o $(LFLAGS)

makenf: makenf.o $(OBJS)
	$(CC) -g -o makenf makenf.o $(OBJS) $(LFLAGS)

makenf_modp: makenf_modp.o $(OBJS)
	$(CC) -g -o makenf_modp makenf_modp.o $(OBJS) $(LFLAGS)

makenf_loop.o: makenf.cc $(OBJS)
	$(CC) -DLOOPER $(CFLAGS) makenf.cc -o makenf_loop.o

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

nflist_loop.o: nflist.cc $(OBJS)
	$(CC) -DLOOPER $(CFLAGS) nflist.cc -o nflist_loop.o

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

looptest: looptest.o looper.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o ratquads.o
	$(CC) -o looptest looptest.o looper.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o ratquads.o $(LFLAGS)

tratquad: tratquad.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o ratquads.o cusp.o mat22.o
	$(CC) -o tratquad tratquad.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o ratquads.o cusp.o mat22.o $(LFLAGS)

homtest: homtest.o $(OBJS)
	$(CC) -o homtest homtest.o $(OBJS) $(LFLAGS)

dimtable_modp: dimtable_modp.o $(OBJS)
	$(CC) -o dimtable_modp dimtable_modp.o $(OBJS) $(LFLAGS)

dimtable: dimtable.o $(OBJS)
	$(CC) -o dimtable dimtable.o $(OBJS) $(LFLAGS)

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

qidltest: qidltest.o primes.o qideal.o qidloop.o quads.o intprocs.o euclid.o geometry.o mat22.o ratquads.o
	$(CC) -o qidltest qidltest.o qidloop.o primes.o qideal.o quads.o intprocs.o euclid.o geometry.o mat22.o  ratquads.o $(LFLAGS)

# DEPENDENCIES
#
# recreate with
# for f in *.cc; do g++ -MM -std=c++11 ${f}; done > t
# and insert t here
#
cusp.o: cusp.cc cusp.h mat22.h ratquads.h quads.h intprocs.h primes.h \
 qideal.h
dimtabeis.o: dimtabeis.cc qidloop.h qideal.h quads.h intprocs.h \
 homspace.h cusp.h mat22.h ratquads.h primes.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
dimtable.o: dimtable.cc qidloop.h qideal.h quads.h intprocs.h homspace.h \
 cusp.h mat22.h ratquads.h primes.h face_relations.h edge_relations.h \
 geometry.h P1N.h hecke.h
dimtable_modp.o: dimtable_modp.cc qidloop.h qideal.h quads.h intprocs.h \
 homspace.h cusp.h mat22.h ratquads.h primes.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
dimtabnew.o: dimtabnew.cc qidloop.h qideal.h quads.h intprocs.h \
 homspace.h cusp.h mat22.h ratquads.h primes.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
dimtabtwist.o: dimtabtwist.cc matprocs.h intprocs.h qidloop.h qideal.h \
 quads.h homspace.h cusp.h mat22.h ratquads.h primes.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h oldforms.h
edge_relations.o: edge_relations.cc mat22.h ratquads.h quads.h intprocs.h \
 primes.h qideal.h homspace.h cusp.h face_relations.h edge_relations.h \
 geometry.h P1N.h hecke.h
euclid.o: euclid.cc euclid.h quads.h intprocs.h geometry.h mat22.h \
 ratquads.h primes.h qideal.h
face_relations.o: face_relations.cc mat22.h ratquads.h quads.h intprocs.h \
 primes.h qideal.h homspace.h cusp.h face_relations.h edge_relations.h \
 geometry.h P1N.h hecke.h
fieldinfo.o: fieldinfo.cc primes.h qideal.h quads.h intprocs.h
geometry.o: geometry.cc geometry.h mat22.h ratquads.h quads.h intprocs.h \
 primes.h qideal.h
hecke.o: hecke.cc hecke.h mat22.h ratquads.h quads.h intprocs.h primes.h \
 qideal.h P1N.h
hecketest.o: hecketest.cc matprocs.h intprocs.h qidloop.h qideal.h \
 quads.h newforms.h ratquads.h oldforms.h primes.h homspace.h cusp.h \
 mat22.h face_relations.h edge_relations.h geometry.h P1N.h hecke.h
hecketest_modp.o: hecketest_modp.cc matprocs.h intprocs.h qidloop.h \
 qideal.h quads.h homspace.h cusp.h mat22.h ratquads.h primes.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
homspace.o: homspace.cc euclid.h quads.h intprocs.h cusp.h mat22.h \
 ratquads.h primes.h qideal.h homspace.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
homtest.o: homtest.cc qidloop.h qideal.h quads.h intprocs.h homspace.h \
 cusp.h mat22.h ratquads.h primes.h face_relations.h edge_relations.h \
 geometry.h P1N.h hecke.h
intprocs.o: intprocs.cc intprocs.h
lf1.o: lf1.cc lf1.h newforms.h ratquads.h quads.h intprocs.h oldforms.h \
 primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
looper.o: looper.cc looper.h quads.h intprocs.h
looptest.o: looptest.cc looper.h quads.h intprocs.h
makenf.o: makenf.cc qidloop.h qideal.h quads.h intprocs.h newforms.h \
 ratquads.h oldforms.h primes.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
makenf_modp.o: makenf_modp.cc qidloop.h qideal.h quads.h intprocs.h \
 newforms.h ratquads.h oldforms.h primes.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
mat22.o: mat22.cc primes.h qideal.h quads.h intprocs.h mat22.h ratquads.h
matprocs.o: matprocs.cc matprocs.h intprocs.h qidloop.h qideal.h quads.h
modularity.o: modularity.cc newforms.h ratquads.h quads.h intprocs.h \
 oldforms.h primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
modularity_modp.o: modularity_modp.cc newforms.h ratquads.h quads.h \
 intprocs.h oldforms.h primes.h qideal.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
moreap1.o: moreap1.cc newforms.h ratquads.h quads.h intprocs.h oldforms.h \
 primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
moreap.o: moreap.cc newforms.h ratquads.h quads.h intprocs.h oldforms.h \
 primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
moreap_loop.o: moreap_loop.cc newforms.h ratquads.h quads.h intprocs.h \
 oldforms.h primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h looper.h
newforms.o: newforms.cc looper.h quads.h intprocs.h newforms.h ratquads.h \
 oldforms.h primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
nflist.o: nflist.cc qidloop.h qideal.h quads.h intprocs.h newforms.h \
 ratquads.h oldforms.h primes.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
nflist_modp.o: nflist_modp.cc qidloop.h qideal.h quads.h intprocs.h \
 newforms.h ratquads.h oldforms.h primes.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
nftest.o: nftest.cc qidloop.h qideal.h quads.h intprocs.h newforms.h \
 ratquads.h oldforms.h primes.h homspace.h cusp.h mat22.h \
 face_relations.h edge_relations.h geometry.h P1N.h hecke.h
oldforms.o: oldforms.cc newforms.h ratquads.h quads.h intprocs.h \
 oldforms.h primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h
P1N.o: P1N.cc P1N.h mat22.h ratquads.h quads.h intprocs.h primes.h \
 qideal.h
P1Ntest.o: P1Ntest.cc looper.h quads.h intprocs.h qidloop.h qideal.h \
 mat22.h ratquads.h primes.h P1N.h
pmanin.o: pmanin.cc newforms.h ratquads.h quads.h intprocs.h oldforms.h \
 primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h looper.h
primes.o: primes.cc primes.h qideal.h quads.h intprocs.h qidloop.h
qideal.o: qideal.cc intprocs.h mat22.h ratquads.h quads.h primes.h \
 qideal.h
qidloop.o: qidloop.cc qidloop.h qideal.h quads.h intprocs.h
qidltest.o: qidltest.cc qidloop.h qideal.h quads.h intprocs.h mat22.h \
 ratquads.h primes.h
quads.o: quads.cc intprocs.h quads.h primes.h qideal.h geometry.h mat22.h \
 ratquads.h
ratquads.o: ratquads.cc ratquads.h quads.h intprocs.h qideal.h mat22.h \
 primes.h geometry.h
roundtest.o: roundtest.cc intprocs.h
testlf1.o: testlf1.cc newforms.h ratquads.h quads.h intprocs.h oldforms.h \
 primes.h qideal.h homspace.h cusp.h mat22.h face_relations.h \
 edge_relations.h geometry.h P1N.h hecke.h lf1.h
tquads.o: tquads.cc looper.h quads.h intprocs.h geometry.h mat22.h \
 ratquads.h primes.h qideal.h
tratquad.o: tratquad.cc ratquads.h quads.h intprocs.h primes.h qideal.h

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
