#Makefile for QUADS directory with test programs

METHOD = 2

GCC=g++
CC = $(GCC)
OPTFLAG = -g -O0 -Wall -DMETHOD=$(METHOD) -fPIC -DNTL_ALL -DUSE_PARI_FACTORING# -DECLIB_MULTITHREAD

# NB If used with a multithreaded build of eclib then you MUST define
# ECLIB_MULTITHREAD above and include the BOOST libraries below.

ECLIB_BASE=/usr/local# $(HOME)/eclib-parallel
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

BOOST_ASIO_LIB = -lboost_system-mt
BOOST_CPPFLAGS =  -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_TR1_UNORDERED_MAP=/\*\*/ -DHAVE_STDCXX_0X=/\*\*/ -DHAVE_UNORDERED_MAP=/\*\*/# -pthread -I/usr/include
BOOST_SYSTEM_LIB = -lboost_system
BOOST_THREAD_LIB = -lboost_thread
BOOST_LDFLAGS = -L/usr/lib $(BOOST_SYSTEM_LIB) $(BOOST_THREAD_LIB)

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

#CFLAGS = -c $(OPTFLAG) $(BOOST_CPPFLAGS) -I$(INCDIR)  -DMETHOD=$(METHOD) -DUSE_XSPLIT
#LFLAGS = -lec -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR) $(BOOST_LDFLAGS)
CFLAGS = -c $(OPTFLAG) -I$(INCDIR)  -DMETHOD=$(METHOD) -DUSE_XSPLIT
LFLAGS = -lec -lntl -lstdc++ -L$(LIBDIR) -Wl,-rpath

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs1 ccs2 ccs3
ccs1: quads.cc fieldinfo.cc cusp.cc homspace.cc homtest.cc hecketest.cc lf1.cc looper.cc looptest.cc
ccs2: moddata.cc modtest.cc mquads.cc newforms.cc oldforms.cc
ccs3: symb.cc symbtest.cc testlf1.cc tmanin.cc tmquads.cc tquads.cc tratquad.cc xtmanin.cc dimtable.cc dimtabeis.cc nftest.cc nflist.cc moreap.cc moreap1.cc

headers:cusp.h homspace.h lf1.h looper.h moddata.h mquads.h newforms.h oldforms.h quads.h ratquads.h symb.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = fieldinfo tquads tratquad looptest modtest symbtest homtest hecketest tmanin moreap moreap1 nftest nflist dimtable dimtabeis # tmquads xtmanin testlf1
tests: $(TESTS)

FIELDS9=1 2 3 7 11 19 43 67 163
FIELDS5=1 2 3 7 11
TESTS9 =  tquads tratquad looptest modtest fieldinfo
TESTS5 =  symbtest homtest hecketest tmanin nftest nflist moreap moreap1 dimtable dimtabeis
export NF_DIR:=nftmp
check: $(TESTS9) $(TESTS5)
	 rm -f t
	 rm -rf $(NF_DIR)
	 mkdir $(NF_DIR)
	 for d in 4 8 3 7 11; do \
         mkdir $(NF_DIR)/2.0.$$d.1; \
	 done
	 LDLIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIBDIR)
	 for p in $(TESTS9); do for d in $(FIELDS9); do \
	 echo "running $$p for d=$$d";\
	 ./$$p < testin/$$p.$$d.in > $$p.$$d.testout 2>/dev/null; diff -a $$p.$$d.testout testout/$$p.$$d.out; \
	 done; done
	 for d in $(FIELDS9); do \
	 echo "running looptest (both conjugates) for d=$$d";\
	 ./looptest < testin/looptest.$${d}a.in > looptest.$$d.testout; diff -a looptest.$$d.testout testout/looptest.$${d}a.out; \
	 done
	 for p in $(TESTS5); do for d in $(FIELDS5); do \
	 echo "running $$p for d=$$d";\
	 ./$$p < testin/$$p.$$d.in > $$p.$$d.testout 2>/dev/null; diff -a $$p.$$d.testout testout/$$p.$$d.out; \
	 done; done
	 rm -rf $(NF_DIR)
	 rm -f *.testout

clean:
	rm -f $(TESTS)
	rm -f *.o *~ *.testout

tquads: tquads.o quads.o
	$(CC) -o tquads tquads.o quads.o $(LFLAGS)

fieldinfo: fieldinfo.o quads.o
	$(CC) -o fieldinfo fieldinfo.o quads.o $(LFLAGS)

tmquads: tmquads.o mquads.o
	$(CC)  -o tmquads tmquads.o mquads.o $(LFLAGS)

OBJS = symb.o moddata.o quads.o looper.o cusp.o homspace.o \
       newforms.o oldforms.o


tmanin: tmanin.o $(OBJS)
	$(CC) -o tmanin tmanin.o $(OBJS) $(LFLAGS)

xtmanin: xtmanin.o $(OBJS)
	$(CC) -o xtmanin xtmanin.o $(OBJS) $(LFLAGS)

testlf1: testlf1.o lf1.o $(OBJS)
	$(CC) -o testlf1 testlf1.o lf1.o $(OBJS) $(LFLAGS)

nftest: nftest.o $(OBJS)
	$(CC) -o nftest nftest.o $(OBJS) $(LFLAGS)

nflist: nflist.o $(OBJS)
	$(CC) -o nflist nflist.o $(OBJS) $(LFLAGS)

moreap: moreap.o $(OBJS)
	$(CC) -o moreap moreap.o $(OBJS) $(LFLAGS)

moreap1: moreap1.o $(OBJS)
	$(CC) -o moreap1 moreap1.o $(OBJS) $(LFLAGS)

looptest: looptest.o looper.o quads.o
	$(CC) -o looptest looptest.o looper.o quads.o $(LFLAGS)

tratquad: tratquad.o quads.o
	$(CC) -o tratquad tratquad.o quads.o $(LFLAGS)

modtest: modtest.o  moddata.o quads.o looper.o
	$(CC) -o modtest modtest.o moddata.o quads.o looper.o $(LFLAGS)

symbtest: symbtest.o symb.o moddata.o quads.o looper.o
	$(CC) -o symbtest symbtest.o symb.o moddata.o quads.o looper.o $(LFLAGS)

homtest: homtest.o symb.o moddata.o quads.o looper.o cusp.o homspace.o
	$(CC) -o homtest homtest.o symb.o moddata.o quads.o looper.o \
                       cusp.o homspace.o $(LFLAGS)

dimtable: dimtable.o symb.o moddata.o quads.o looper.o cusp.o homspace.o
	$(CC) -o dimtable dimtable.o symb.o moddata.o quads.o looper.o \
                       cusp.o homspace.o $(LFLAGS)

dimtabeis: dimtabeis.o symb.o moddata.o quads.o looper.o cusp.o homspace.o
	$(CC) -o dimtabeis dimtabeis.o symb.o moddata.o quads.o looper.o \
                       cusp.o homspace.o $(LFLAGS)

hecketest: hecketest.o symb.o moddata.o quads.o looper.o cusp.o homspace.o
	$(CC) -o hecketest hecketest.o symb.o moddata.o quads.o looper.o \
                       cusp.o homspace.o $(LFLAGS)

roundtest: roundtest.o quads.o
	$(CC) -o roundtest roundtest.o quads.o $(LFLAGS)

# DEPENDENCIES

quads.o: quads.cc quads.h
fieldinfo.o: quads.h
ratquads.o: quads.h ratquads.h
tratquad.o: tratquad.cc ratquads.h quads.h
tquads.o: tquads.cc quads.h
mquads.o: mquads.cc mquads.h
tmquads.o: tmquads.cc mquads.h
modtest.o: modtest.cc moddata.h quads.h looper.h
symbtest.o: symbtest.cc symb.h moddata.h ratquads.h quads.h looper.h
homtest.o: homtest.cc cusp.h homspace.h symb.h moddata.h ratquads.h quads.h \
           looper.h
dimtable.o: dimtable.cc cusp.h homspace.h symb.h moddata.h ratquads.h quads.h \
           looper.h
dimtabeis.o: dimtabeis.cc cusp.h homspace.h symb.h moddata.h ratquads.h \
	quads.h looper.h
hecketest.o: hecketest.cc cusp.h homspace.h symb.h moddata.h ratquads.h \
             quads.h looper.h
tmanin.o: tmanin.cc looper.h symb.h moddata.h ratquads.h quads.h \
          looper.h homspace.h newforms.h oldforms.h cusp.h
testlf1.o : testlf1.cc oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h lf1.h
nftest.o : nftest.cc oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h
nflist.o : nflist.cc oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h
moreap.o : moreap.cc oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h
moreap1.o : moreap1.cc oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h
moddata.o: moddata.cc moddata.h quads.h
symb.o: symb.cc symb.h moddata.h ratquads.h quads.h
cusp.o: quads.h ratquads.h moddata.h symb.h cusp.h cusp.cc
homspace.o: quads.h ratquads.h moddata.h symb.h cusp.h homspace.h homspace.cc
oldforms.o: oldforms.h oldforms.cc moddata.h quads.h
newforms.o: newforms.h newforms.cc moddata.h quads.h symb.h homspace.h
looper.o: looper.h looper.cc quads.h
looptest.o: looper.h looptest.cc quads.h
lf1.o: lf1.cc lf1.h moddata.h quads.h symb.h homspace.h ratquads.h
test.o: test.cc


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
