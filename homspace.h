// FILE HOMSPACE.H: Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include "arith.h"
#include "method.h"
#include "subspace.h"
#include "smatrix.h"
#ifdef USE_SMATS
#include "smatrix_elim.h"
#endif
#include "moddata.h"
#include "cusp.h"
#include "symb.h"

class homspace :public symbdata {
friend class newforms;
private:
  int verbose;
   int *coordindex,*needed,*freegens;
   long rk,denom1,denom2;
   subspace kern;
   modsym *freemods;
   mat opmat(long, int dual=1, int verb=0);
   mat opmat_restricted(int i,const subspace& s, int dual, int verb=0);
  // versions returning an smat:
   smat s_opmat(int i,int dual,int verb=0);
   smat s_opmat_restricted(int i,const ssubspace& s, int dual,int verb=0);
   long matdim(void) {return dimension;} 
   long matden(void) {return denom3;}
   vector<long> eigrange(long i);
 protected:
   mat coord, projcoord;
   long dimension, denom3, ncusps;
   mat relmat; long numrel, maxnumrel;
   void userel(vec& rel);
 public:
   homspace(const Quad& n, int hp, int verb);
   ~homspace() {delete[] coordindex; delete[] needed; 
                delete[] freegens; delete[] freemods;
              }
   long h1dim() const {return dimension;}  // No confusion with subspace::dim
   long h1denom() const {return denom3;}
   long h1ncusps() const {return ncusps;}
 public:
   vec chain(const symb& s) const;
   vec chaincd(const Quad& c, const Quad& d) const;
   vec projchaincd(const Quad& c, const Quad& d) const;
   vec chain(const Quad& nn, const Quad& dd) const;
   vec chain(const RatQuad& r) const 
        {return chain(num(r),den(r));}
   vec kernelpart(const vec& v) const 
        {return v[pivots(kern)];}
   vec cycle(const symb& s) const 
        {return kernelpart(chain(s));}
   vec cycle(const Quad& n, const Quad& d) const 
        {return kernelpart(chain(n,d));}
   vec cycle(const RatQuad& r) const 
        {return kernelpart(chain(num(r),den(r)));}
   vec cycle(const modsym& m) const
     {return cycle(m.beta())-cycle(m.alpha());}
   vec projcycle(const Quad& n, const Quad& d) const;
   vec applyop(const matop& mlist, const RatQuad& q) const;
   vec applyop(const matop& mlist, const modsym& m) const 
                   {return applyop(mlist,m.beta())-applyop(mlist,m.alpha());} 
   mat calcop(const string opname, const Quad& p, const matop& mlist, int dual=1, int display=0) const;
   smat s_calcop(const string  opname, const Quad& p, const matop& mlist, 
		 int dual, int display) const;
   mat calcop_restricted(const string opname, const Quad& p, const matop& mlist, const subspace& s, int dual, int display) const;
   smat s_calcop_restricted(const string opname, const Quad& p, const matop& mlist, const ssubspace& s, int dual, int display) const;
public:
   mat heckeop(const Quad& p, int dual=1, int display=0) const;
   smat s_heckeop(const Quad& p, int dual, int display) const;
   mat heckeop_restricted(const Quad& p, const subspace& s, int dual, int display) const;
   smat s_heckeop_restricted(const Quad& p, const ssubspace& s, int dual, int display) const;
   mat wop(const Quad& q, int dual=1, int display=0) const;
   mat fricke(int dual=1, int display=0) const;
//   mat conj(int display=0) const;
   vec maninvector(const Quad& p) const;
   vec projmaninvector(const Quad& p) const;
   vec manintwist(const Quad& lambda, const vector<Quad>& res, int* chitable) const;
   vec newhecke(const Quad& p, const Quad& n, const Quad& d) const;
};

#ifndef USE_XSPLIT

int startswith(const vector<long>& a, const vector<long>& b, long l);

class splitter {
private:
  homspace* h1;
  long maxdepth, depth, h1denom;
  subspace* nest;
  vector<long> aplist; Quadlist plist;
  vec basisvector;
  int *havemat;
  mat *opmats;
  int dual, method, verbose;
//method=0 for usual elimination using longs
//method=1 uses long longs
//method=2 uses large prime modulus P for elimination only
//method=3 works mod P throughout, lifting only at end.
public:
  splitter(homspace* h, long d, int dualflag, int meth, int v=0);
  ~splitter(void) {delete[] nest; delete[] havemat; delete[] opmats;}
  void splitoff(const vector<long>& apl);
  vec getbasis() const {return basisvector;}
};

#endif  // #ifndef USE_XSPLIT

#endif
