// FILE MODDATA.H: Declaration of class moddata  QUADS version

#if     !defined(_MODDATA_H)
#define _MODDATA_H      1       //flags that this file has been included

#include "quads.h"

class level {
public:
  Quad modulus, conj_modulus;
  int plusflag;
  vector<Quad> plist,dlist,primelist;
  long npdivs,ndivs,normod,nap;
  int is_square;
  int is_Galois_stable;

protected:
//If modulus=(a,b) with norm normod, n0=gcd(a,b), n0m0=normod/n0=n0*m0
  long n0m0, n0,m0, wmodz;
  Quad reduce(const Quad& a) const;  // return a reduced mod modulus
  long numres(const Quad& a) const;  // what number is this residue a mod modulus?
 Quad resnum(long i) const;  // which is the i'th residue mod modulus?
// The constructor:
  explicit level(const Quad& n, long neigs=20); // neigs controls the maximum depth for recursion in newform finding
};

class moddata :public level {
protected:
 long phi,psi,nsymb,nsymb1,nsymb2;
 vector<long> invlist,dstarts;
 vector<long> noninvlist,noninvdlist;
 long code(const Quad& res) const {return invlist[numres(res)];}
public:
 explicit moddata(const Quad& n);                                //constructor
 ~moddata();                                    //destructor
 void display() const;
 int check(int verbose=0) const;     //checks whether resnum & numres work OK
 void abort(char* mess) const {cerr<<"Out of memory ("<<mess<<")\n";  exit(1);}
};

#endif
