#Makefile for JQUADS directory with test programs

MAXCLASSNUM = 2
METHOD = 2

# Include machine-dependent directories
include make.inc.fermat

#########################################################################
# The next line is needed so that suffix rules using .cc can be defined #
#########################################################################
.SUFFIXES: .cc $(SUFFIXES)

CC = gcc

CFLAGS = -g -c -O0  -I$(INCDIR) -I$(NTLINCDIR) -DMAX_CLASSNUM=$(MAXCLASSNUM) -DMETHOD=$(METHOD) -DUSE_XSPLIT -Dtestbezout -DNTL_ALL #-DUSE_SMATS
LFLAGS =  -L$(LIBDIR) -ljcntl -L$(NTLLIBDIR) -lntl -L$(HOME)/sage-current/local/lib -lpari -lgmp -lstdc++

######################
# GENERAL MAKE RULES #
######################

# COMPILE .cc to .o but do not link
.cc.o:
	@echo Using my own suffix rule .cc.o for compile step...
	$(CC) $(CFLAGS) $<
#	$(CC) $(CFLAGS) -c $< -o $@

# LINK ONLY (how to make an executable from a .o)
.o:
	@echo Using my own .o suffix rule for link step...
	$(CC) -o $@ $< $(LFLAGS)

# COMPILE AND LINK
#%: %.cc
#	$(CC) $(CFLAGS) $(LFLAGS) $< -o $@

#########################
# INDIVIDUAL MAKE RULES #
#########################

sources: ccs headers
	chmod a+r *.h *.cc

ccs: ccs1 ccs2 ccs3
ccs1: homspace.cc homtest.cc hecketest.cc lf1.cc looper.cc looptest.cc
ccs2: manin.cc moddata.cc modtest.cc mquads.cc newforms.cc oldforms.cc 
ccs3: quads.cc symb.cc symbtest.cc testlf1.cc tmanin.cc tmquads.cc tquads.cc xtmanin.cc

headers: homspace.h lf1.h looper.h manin.h moddata.h mquads.h newforms.h oldforms.h quads.h symb.h

TESTS = tquads tmquads looptest qidltest modtest symbtest homtest hecketest tmanin xtmanin testlf1 nftest

JTESTS = tquads looptest qidltest modtest symbtest homtest hecketest

tests: $(TESTS)
jtests: $(JTESTS)

jpublic: jtests
	chmod a+rx $(JTESTS)
	chmod a+r *.h *.cc Makefile

shar: sources
	rm -f sharfile
	shar -i files -o sharfile
	compress sharfile

notests:
	rm -f $(TESTS)

clean: notests
	rm -f *.o

cleaner: clean
	rm -f *~
#	sccs clean

tquads: quads.o tquads.o
	$(CC) -o tquads quads.o tquads.o $(LFLAGS)

tmquads: tmquads.o mquads.o
	$(CC)  -o tmquads tmquads.o mquads.o $(LFLAGS)

looptest: looptest.o looper.o quads.o
	$(CC) -o looptest looptest.o looper.o quads.o $(LFLAGS)

qidltest: qidltest.o qidloop.o quads.o
	$(CC) -o qidltest qidltest.o qidloop.o quads.o $(LFLAGS)

modtest: modtest.o moddata.o qidloop.o looper.o quads.o
	$(CC) -o modtest modtest.o moddata.o quads.o qidloop.o looper.o $(LFLAGS)

symbtest: symbtest.o symb.o moddata.o qidloop.o looper.o quads.o
	$(CC) -o symbtest symbtest.o symb.o moddata.o quads.o qidloop.o looper.o $(LFLAGS)

homtest: homtest.o homspace.o symb.o moddata.o qidloop.o looper.o quads.o
	$(CC) -o homtest homtest.o homspace.o symb.o moddata.o quads.o qidloop.o looper.o $(LFLAGS)

hecketest: hecketest.o homspace.o symb.o moddata.o qidloop.o looper.o quads.o
	$(CC) -o hecketest hecketest.o homspace.o symb.o moddata.o qidloop.o \
                           looper.o quads.o $(LFLAGS)

roundtest: roundtest.o quads.o
	$(CC) -o roundtest roundtest.o quads.o $(LFLAGS)


OBJS = symb.o moddata.o quads.o looper.o homspace.o \
       newforms.o oldforms.o manin.o

tmanin: tmanin.o $(OBJS)
	$(CC) -o tmanin tmanin.o $(OBJS) $(LFLAGS)

xtmanin: xtmanin.o $(OBJS)
	$(CC) -o xtmanin xtmanin.o $(OBJS) $(LFLAGS)

testlf1: testlf1.o lf1.o $(OBJS)
	$(CC) -o testlf1 testlf1.o lf1.o $(OBJS) $(LFLAGS)

nftest: nftest.o $(OBJS)
	$(CC) -o nftest nftest.o $(OBJS) $(LFLAGS)

# DEPENDENCIES

quadarith.o: quadarith.cc quadarith.h field.h
field.o: field.cc field.h quadarith.h
qideal.o: qideal.cc qideal.h quadarith.h field.h
primes.o: primes.cc primes.h primes1.cc primes2.cc qideal.h quadarith.h field.h
quads.o: quads.cc quads.h field.cc field.h quadarith.cc quadarith.h qideal.cc qideal.h primes.cc primes.h primes1.cc primes2.cc
quads.h: field.h quadarith.h qideal.h primes.h

tquads.o: tquads.cc quads.h
mquads.o: mquads.cc mquads.h
tmquads.o: tmquads.cc mquads.h
modtest.o: modtest.cc moddata.h quads.h qidloop.h
symbtest.o: symbtest.cc symb.h moddata.h quads.h qidloop.h
homtest.o: homtest.cc homspace.h symb.h moddata.h quads.h qidloop.h looper.h
hecketest.o: hecketest.cc homspace.h symb.h moddata.h qidloop.h quads.h
tmanin.o: tmanin.cc manin.h looper.h symb.h moddata.h quads.h \
          looper.h homspace.h newforms.h oldforms.h
testlf1.o : testlf1.cc manin.h oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h symb.h lf1.h 
nftest.o : nftest.cc manin.h oldforms.h moddata.h quads.h newforms.h \
            looper.h homspace.h symb.h 
moddata.o: moddata.cc looper.h moddata.h quads.h
symb.o: symb.cc symb.h moddata.h quads.h
homspace.o: quads.h moddata.h symb.h homspace.h homspace.cc
manin.o: manin.h manin.cc quads.h moddata.h symb.h \
         homspace.h looper.h newforms.h
oldforms.o: oldforms.h oldforms.cc moddata.h quads.h
newforms.o: newforms.h newforms.cc moddata.h quads.h symb.h homspace.h
looper.o: looper.h looper.cc quads.h
looptest.o: looper.h looptest.cc quads.h
qidloop.o: qidloop.h qidloop.cc quads.h
qidltest.o: qidloop.h qidltest.cc quads.h
lf1.o: lf1.cc lf1.h moddata.h quads.h symb.h homspace.h
test.o: test.cc


