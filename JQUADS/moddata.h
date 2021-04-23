// FILE moddata.h: Declaration of class moddata  Q(sqrt(-5)) version

#if     !defined(_MODDATA_H)
#define _MODDATA_H      1       //flags that this file has been included

#include "arith.h"   // even better, get quads.h to include this
#include "quads.h"
#include "qideal.h"

class set_of_residues {
  long c,ac,bc;
  Quadlist resl;
  Quad box_resnum(long i) const;   // which is the i'th residue mod modulus?
protected:
  long normod;
public:
  set_of_residues(const Qideal& n);   // constructor

  // member functions
  long numres(const Quad& a) const;  // what number is residue a mod modulus?
  Quad resnum(long i) const;         // which is the i'th residue mod modulus?
  int check(int verbose=0) const;    // checks whether resnum & numres work OK
  void display(ostream&s = cout, int verbose=0) const;
};

class quotient_ring : public set_of_residues {
  vector<long> invlist;
public:
  quotient_ring(const Qideal& n, long phi_n);          //constructor
  void display(ostream& s = cout, int verbose=0) const;
protected:
  vector<long> noninvlist;
  long code(const Quad& res) const {return invlist[numres(res)];}
//  Quad noninvresnum(long i) const {return resnum(noninvlist[i]);}
};

class moddata {
// problematic stuff which used to be in the ancestor "level"
//    public:
//    static Quadlist plist,dlist,primelist;
//    static long npdivs,ndivs,normod,nap;
//    level(const Quad& n, long neigs=10);              //constructor
public:
  Qideal modulus;
  prime_factn nfactn;
protected:
  long phi,psi;
public:
  moddata(const Qideal& n);                                //constructor
  void display(ostream& s = cout) const;
  long get_phi() const {return phi;}
};

// some other general-purpose functions (could be in quads itself)

// #include <values.h>
// inline double psif(complex z) {  return cos(4*PI*real(z));}
// inline double psig(complex z) {  return sin(4*PI*real(z));}
// int squaremod(const Quad& a, const Quad& m, const Quadlist& reslist);
// int* makechitable(const Quad& lambda, const Quadlist& reslist);
// double gauss(const Quad& m, const Quadlist& reslist);

#endif

// END OF FILE moddata.h
