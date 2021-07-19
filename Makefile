#Makefile for bianchi-progs with test programs

# eclib is a requirement.  If installed in the usual place /usr/local
# the following line might not be necessary.  If installed anywhere
# else set the install directory here:

#ECLIB_BASE=/usr/local
ECLIB_BASE=$(HOME)/eclib
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

USE_SMATS=1
#  It is recommended to USE_SMATS when finding rational newforms
#  only, since the linear algebra is done correctly modulo p=2^30-35
#  and only lifted to Z when an eigenvector is found. Programs which
#  compute Hecke matrices and their characteristic polynomials such
#  as hecketest.cc try to lift the whole modular symbol space and
#  this easily fails.  In that case do not USE_SMATS and the linear
#  algebra will be done using multiprecision integer arithmetic
#  instead.  This is much slower: e.g. with level (128), field 1 it
#  takes 20m instead of <1s.

GCC=g++ -std=c++11
CC = $(GCC)
OPTFLAG = -O3 -Wall -fPIC

# NB If used with a multithreaded build of eclib then you MUST define
# USE_BOOST=1 below so that the correct compiler and linker stuff is
# appended below.  Otherwise set USE_BOOST=0.

ifeq ($(USE_BOOST), 1)
 BOOST_ASIO_LIB = -lboost_system-mt
 BOOST_CPPFLAGS =   -DECLIB_MULTITHREAD -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_TR1_UNORDERED_MAP=/\*\*/ -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_UNORDERED_MAP=/\*\*/# -pthread -I/usr/include
 BOOST_SYSTEM_LIB = -lboost_system
 BOOST_THREAD_LIB = -lboost_thread
 BOOST_LDFLAGS = -L/usr/lib $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)
endif

# for profiling:
#CFLAGS = -c -pg $(OPTFLAG) $(BOOST_CPPFLAGS) -I$(INCDIR) -DUSE_SMATS=$(USE_SMATS)
#LFLAGS = -pg -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

#for normal use:
CFLAGS = -c -g $(OPTFLAG) $(BOOST_CPPFLAGS) -I$(INCDIR) -DUSE_SMATS=$(USE_SMATS)
LFLAGS = -lec -lntl -lstdc++  -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)

all: tests

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs1 ccs2 ccs3 ccs4
ccs1: intprocs.cc quads.cc mat22.cc fieldinfo.cc cusp.cc homtest.cc hecketest.cc lf1.cc looper.cc looptest.cc euclid.cc geometry.cc
ccs2: moddata.cc modtest.cc mquads.cc newforms.cc oldforms.cc homspace.cc edge_relations.cc face_relations.cc hecke.cc
ccs3: symb.cc symbtest.cc testlf1.cc tmanin.cc pmanin.cc tmquads.cc tquads.cc tratquad.cc xtmanin.cc dimtable.cc dimtabeis.cc nftest.cc nflist.cc moreap.cc moreap1.cc moreap_loop.cc modularity.cc modularity_modp.cc
ccs4: qideal.cc qidloop.cc primes.cc qidltest.cc

headers: intprocs.h cusp.h homspace.h lf1.h looper.h moddata.h mquads.h newforms.h oldforms.h quads.h ratquads.h symb.h euclid.h geometry.h qideal.h primes.h qidloop.h mat22.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = fieldinfo tquads qidltest tratquad looptest modtest symbtest homtest hecketest tmanin moreap moreap1 nftest nflist dimtable dimtabeis modularity modularity_modp # tmquads xtmanin testlf1
tests: $(TESTS)

DISCS9=4 8  3 7 11 19 43 67 163
FIELDS9=1 2 3 7 11 19 43 67 163
FIELDS7=1 2 3 7 11 19 43
FIELDS5=1 2 3 7 11
FIELDS1=1
TESTS9 =  tquads tratquad looptest modtest fieldinfo symbtest homtest dimtable dimtabeis hecketest tmanin nftest nflist moreap moreap1 modularity modularity_modp
ALL_TESTS = $(TESTS9)

test_input_dir = testin
test_output_dir = testout

check_run = echo -n "Testing $${prog} for d=$${d}..."; ./$${prog} < $(test_input_dir)/$${prog}.$${d}.in > $${prog}.$${d}.testout 2>/dev/null && if diff -q $${prog}.$${d}.testout $(test_output_dir)/$${prog}.$${d}.out; then echo "$${prog} for d=$${d} completed successfully"; else echo "$${prog} for d=$${d} failed"; diff $${prog}.$${d}.testout $(test_output_dir)/$${prog}.$${d}.out; fi || exit $$?


export NF_DIR:=nftmp
check: $(ALL_TESTS)
	 @echo Test setup: create temporary directories
	 rm -f t
	 rm -rf $(NF_DIR)
	 mkdir $(NF_DIR)
	 for d in $(DISCS9); do mkdir $(NF_DIR)/2.0.$$d.1; done
	 @echo
	 @echo running $(TESTS9) on fields $(FIELDS9)...
	 @for prog in $(TESTS9); do for d in $(FIELDS9); do $(check_run); done; done
	 @echo
	 @echo running looptest \(both conjugates\) on fields $(FIELDS9)...
	 @for d in $(FIELDS9); do for prog in looptest; do $(check_run); done; done
	 @echo
	 @echo Tidy up: remove temporary directories and output test files
	 rm -rf $(NF_DIR)
	 rm -f *.testout
	 @echo Tests completed

clean:
	rm -f $(TESTS)
	rm -f *.o *~ *.testout

tquads: tquads.o quads.o looper.o intprocs.o euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o
	$(CC) -o tquads tquads.o quads.o looper.o intprocs.o euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o $(LFLAGS)

fieldinfo: fieldinfo.o quads.o euclid.o geometry.o intprocs.o qideal.o qidloop.o primes.o mat22.o
	$(CC) -o fieldinfo fieldinfo.o quads.o euclid.o geometry.o intprocs.o qideal.o qidloop.o primes.o mat22.o $(LFLAGS)

tmquads: tmquads.o mquads.o
	$(CC)  -o tmquads tmquads.o mquads.o $(LFLAGS)

OBJS = symb.o moddata.o quads.o intprocs.o euclid.o geometry.o looper.o cusp.o homspace.o \
       newforms.o oldforms.o edge_relations.o face_relations.o hecke.o qideal.o qidloop.o \
       primes.o mat22.o ratquads.o

tmanin: tmanin.o $(OBJS)
	$(CC) -g -o tmanin tmanin.o $(OBJS) $(LFLAGS)

tmanin_loop.o: tmanin.cc $(OBJS)
	$(CC) -DLOOPER $(CFLAGS) tmanin.cc -o tmanin_loop.o

tmanin_loop: tmanin_loop.o $(OBJS)
	$(CC) -g -o tmanin_loop tmanin_loop.o $(OBJS) $(LFLAGS)

pmanin: pmanin.o $(OBJS)
	$(CC) -o pmanin pmanin.o $(OBJS) $(LFLAGS)

xtmanin: xtmanin.o $(OBJS)
	$(CC) -o xtmanin xtmanin.o $(OBJS) $(LFLAGS)

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

looptest: looptest.o looper.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o
	$(CC) -o looptest looptest.o looper.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o mat22.o $(LFLAGS)

tratquad: tratquad.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o ratquads.o
	$(CC) -o tratquad tratquad.o quads.o intprocs.o  euclid.o geometry.o qideal.o qidloop.o primes.o ratquads.o $(LFLAGS)

modtest: modtest.o $(OBJS)
	$(CC) -o modtest modtest.o $(OBJS) $(LFLAGS)

symbtest: symbtest.o $(OBJS)
	$(CC) -o symbtest symbtest.o $(OBJS) $(LFLAGS)

homtest: homtest.o $(OBJS)
	$(CC) -o homtest homtest.o $(OBJS) $(LFLAGS)

dimtable: dimtable.o $(OBJS)
	$(CC) -o dimtable dimtable.o $(OBJS) $(LFLAGS)

dimtabeis: dimtabeis.o $(OBJS)
	$(CC) -o dimtabeis dimtabeis.o $(OBJS) $(LFLAGS)

hecketest: hecketest.o $(OBJS)
	$(CC) -o hecketest hecketest.o $(OBJS) $(LFLAGS)

roundtest: roundtest.o quads.o
	$(CC) -o roundtest roundtest.o quads.o $(LFLAGS)

qidltest: qidltest.o primes.o qideal.o qidloop.o quads.o intprocs.o euclid.o geometry.o
	$(CC) -o qidltest qidltest.o qidloop.o primes.o qideal.o quads.o intprocs.o euclid.o geometry.o $(LFLAGS)

# DEPENDENCIES
#
# recreate with
# for f in *.cc; do g++ -MM -std=c++11 ${f}; done
#
intprocs.o: intprocs.h intprocs.cc
cusp.o: cusp.cc cusp.h moddata.h quads.h ratquads.h
ratquads.o: ratquads.cc quads.h ratquads.h qideal.h
dimtabeis.o: dimtabeis.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h looper.h
dimtable.o: dimtable.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h looper.h
fieldinfo.o: fieldinfo.cc quads.h intprocs.h qideal.h qidloop.h
hecketest.o: hecketest.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h
homspace.o: homspace.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h euclid.h geometry.h
edge_relations.o: edge_relations.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h geometry.h
face_relations.o: face_relations.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h geometry.h
hecke.o: hecke.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h
homtest.o: homtest.cc homspace.h moddata.h quads.h ratquads.h cusp.h symb.h looper.h
lf1.o: lf1.cc lf1.h newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
looper.o: looper.cc looper.h quads.h
looptest.o: looptest.cc looper.h quads.h
moddata.o: moddata.cc moddata.h quads.h ratquads.h
modtest.o: modtest.cc moddata.h quads.h ratquads.h looper.h
modularity.o: modularity.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
modularity_modp.o: modularity_modp.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
moreap1.o: moreap1.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
moreap.o: moreap.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
moreap_loop.o: moreap_loop.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h looper.h
mquads.o: mquads.cc mquads.h
newforms.o: newforms.cc looper.h quads.h oldforms.h moddata.h ratquads.h newforms.h homspace.h cusp.h symb.h
nflist.o: nflist.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
nftest.o: nftest.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
oldforms.o: oldforms.cc oldforms.h moddata.h quads.h ratquads.h newforms.h homspace.h cusp.h symb.h
pmanin.o: pmanin.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h looper.h
quads.o: quads.cc quads.h intprocs.h geometry.h
euclid.o: euclid.cc euclid.h quads.h geometry.h mat22.h
geometry.o: geometry.cc geometry.h quads.h mat22.h
symb.o: symb.cc symb.h moddata.h quads.h ratquads.h euclid.h geometry.h
symbtest.o: symbtest.cc symb.h moddata.h quads.h ratquads.h looper.h
testlf1.o: testlf1.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h lf1.h
tmanin.o: tmanin.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h
tmquads.o: tmquads.cc mquads.h
tquads.o: tquads.cc quads.h looper.h geometry.h intprocs.h qideal.h qidloop.h primes.h
tratquad.o: tratquad.cc ratquads.cc ratquads.h quads.h
xtmanin.o: xtmanin.cc newforms.h oldforms.h moddata.h quads.h ratquads.h homspace.h cusp.h symb.h

mat22.o: mat22.cc mat22.h quads.h
qideal.o: qideal.cc qideal.h quads.h mat22.h intprocs.h
primes.o: primes.cc primes.h qideal.h quads.h mat22.h intprocs.h
qidloop.o: qidloop.cc qidloop.h qideal.h quads.h mat22.h
qidltest.o: qidltest.cc qidloop.h primes.h qideal.h quads.h mat22.h

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
