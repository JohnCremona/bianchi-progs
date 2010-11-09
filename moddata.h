// FILE MODDATA.H: Declaration of class moddata  QUADS version

#if     !defined(_MODDATA_H)
#define _MODDATA_H      1       //flags that this file has been included

#include "arith.h"   // even better, get quads.h to include this
#include "quads.h"

class level {
public:
 static Quad modulus;
 static int plusflag;
 static vector<Quad> plist,dlist,primelist;
 static long npdivs,ndivs,normod,nap;
protected:
//If modulus=(a,b) with norm normod, n0=gcd(a,b), n0m0=normod/n0=n0*m0
 static long n0m0, n0,m0, wmodz;
 static long numres(const Quad& a);  // what number is this residue a mod modulus?
 static Quad resnum(long i);  // which is the i'th residue mod modulus?
// The constructor:
 level(const Quad& n, long neigs=10); 
};

class moddata :public level {
protected:
 long phi,psi,nsymb,nsymb1,nsymb2;
 vector<long> invlist,dstarts;
 vector<long> noninvlist,noninvdlist;
 long code(const Quad& res) const {return invlist[numres(res)];}
public:
 moddata(const Quad& n);                                //constructor
 ~moddata();                                    //destructor
 void display() const;
 int check(int verbose=0) const;     //checks whether resnum & numres work OK
 void abort(char* mess) const {cerr<<"Out of memory ("<<mess<<")\n";  exit(1);}
};

#endif
