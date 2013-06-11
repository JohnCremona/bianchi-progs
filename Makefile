#Makefile for QUADS directory with test programs

METHOD = 2

GCC=g++
CC = $(GCC)
OPTFLAG = -g  -Wall -DMETHOD=$(METHOD) -DNTL_ALL -DUSE_PARI_FACTORING

ECLIB_BASE=$(HOME)/eclib
INCDIR = $(ECLIB_BASE)/include
LIBDIR = $(ECLIB_BASE)/lib

CLEAN = rcsclean
RANLIB = ranlib
CP = cp -p

CFLAGS = -c $(OPTFLAG)  -I$(INCDIR)  -DMETHOD=$(METHOD) -DUSE_XSPLIT
LFLAGS = -lec -L$(LIBDIR) -Wl,-rpath -Wl,$(LIBDIR)

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs1 ccs2 ccs3
ccs1: cusp.cc homspace.cc homtest.cc hecketest.cc lf1.cc looper.cc looptest.cc
ccs2: manin.cc moddata.cc modtest.cc mquads.cc newforms.cc oldforms.cc
ccs3: quads.cc symb.cc symbtest.cc testlf1.cc tmanin.cc tmquads.cc tquads.cc tratquad.cc xtmanin.cc dimtable.cc

headers:cusp.h homspace.h lf1.h looper.h manin.h moddata.h mquads.h newforms.h oldforms.h quads.h ratquads.h symb.h

%.o:   %.cc
	$(CC) $(CFLAGS) $<

TESTS = tquads tratquad looptest modtest symbtest homtest hecketest tmanin  nftest # tmquads xtmanin testlf1
tests: $(TESTS)

FIELDS9=1 2 3 7 11 19 43 67 163
FIELDS5=1 2 3 7 11
TESTS9 =  tquads tratquad looptest modtest
TESTS5 =  symbtest homtest hecketest tmanin dimtable
check: $(TESTS9) $(TESTS5)
	 rm -f t
	 LDLIBRARY_PATH=$(LD_LIBRARY_PATH):$(LIBDIR)
	 for p in $(TESTS9); do for d in $(FIELDS9); do \
	 echo "running $$p for d=$$d";\
	 ./$$p < testin/$$p.$$d.in > $$p.$$d.testout && diff -a $$p.$$d.testout testout/$$p.$$d.out; \
	 done; done
	 for d in $(FIELDS9); do \
	 echo "running looptest (both conjugates) for d=$$d";\
	 ./looptest < testin/looptest.$${d}a.in > looptest.$$d.testout && diff -a looptest.$$d.testout testout/looptest.$${d}a.out; \
	 done
	 for p in $(TESTS5); do for d in $(FIELDS5); do \
	 echo "running $$p for d=$$d";\
	 ./$$p < testin/$$p.$$d.in > $$p.$$d.testout && diff -a $$p.$$d.testout testout/$$p.$$d.out; \
	 done; done
	 rm -f *.testout

clean:
	rm -f $(TESTS)
	rm -f *.o *~ *.testout

tquads: tquads.o quads.o
	$(CC) -o tquads tquads.o quads.o $(LFLAGS)

tmquads: tmquads.o mquads.o
	$(CC)  -o tmquads tmquads.o mquads.o $(LFLAGS)

OBJS = symb.o moddata.o quads.o looper.o cusp.o homspace.o \
       newforms.o oldforms.o manin.o


tmanin: tmanin.o $(OBJS)
	$(CC) -o tmanin tmanin.o $(OBJS) $(LFLAGS)

xtmanin: xtmanin.o $(OBJS)
	$(CC) -o xtmanin xtmanin.o $(OBJS) $(LFLAGS)

testlf1: testlf1.o lf1.o $(OBJS)
	$(CC) -o testlf1 testlf1.o lf1.o $(OBJS) $(LFLAGS)

nftest: nftest.o $(OBJS)
	$(CC) -o nftest nftest.o $(OBJS) $(LFLAGS)

looptest: looptest.o looper.o quads.o
	$(CC) -o looptest looptest.o looper.o quads.o $(LFLAGS)

tratquad: tratquad.o
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

hecketest: hecketest.o symb.o moddata.o quads.o looper.o cusp.o homspace.o
	$(CC) -o hecketest hecketest.o symb.o moddata.o quads.o looper.o \
                       cusp.o homspace.o $(LFLAGS)

roundtest: roundtest.o quads.o
	$(CC) -o roundtest roundtest.o quads.o $(LFLAGS)

# DEPENDENCIES

quads.o: quads.cc quads.h
ratquads.o: quads.h ratquads.h
tratquad.o: tratquad.cc ratquads.h quads.h
tquads.o: tquads.cc quads.h
mquads.o: mquads.cc mquads.h
tmquads.o: tmquads.cc mquads.h
modtest.o: modtest.cc moddata.h quads.h looper.h
symbtest.o: symbtest.cc symb.h moddata.h ratquads.h quads.h looper.h
homtest.o: homtest.cc cusp.h homspace.h symb.h moddata.h ratquads.h quads.h \
           looper.h
dimtable.o: homtest.cc cusp.h homspace.h symb.h moddata.h ratquads.h quads.h \
           looper.h
hecketest.o: hecketest.cc cusp.h homspace.h symb.h moddata.h ratquads.h \
             quads.h looper.h
tmanin.o: tmanin.cc manin.h looper.h symb.h moddata.h ratquads.h quads.h \
          looper.h homspace.h newforms.h oldforms.h cusp.h
testlf1.o : testlf1.cc manin.h oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h lf1.h
nftest.o : nftest.cc manin.h oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h cusp.h ratquads.h symb.h
moddata.o: moddata.cc moddata.h quads.h
symb.o: symb.cc symb.h moddata.h ratquads.h quads.h
cusp.o: quads.h ratquads.h moddata.h symb.h cusp.h cusp.cc
homspace.o: quads.h ratquads.h moddata.h symb.h cusp.h homspace.h homspace.cc
manin.o: manin.h manin.cc quads.h ratquads.h moddata.h symb.h cusp.h \
         homspace.h looper.h newforms.h
oldforms.o: oldforms.h oldforms.cc moddata.h quads.h
newforms.o: newforms.h newforms.cc moddata.h quads.h symb.h homspace.h
looper.o: looper.h looper.cc quads.h
looptest.o: looper.h looptest.cc quads.h
lf1.o: lf1.cc lf1.h moddata.h quads.h symb.h homspace.h ratquads.h
test.o: test.cc


