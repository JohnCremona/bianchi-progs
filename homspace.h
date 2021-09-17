// FILE HOMSPACE.H: Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include "cusp.h"
#include "face_relations.h"

class homspace {
friend class newforms;
public:
  int verbose;
  int cuspidal;  // if 1 then compute cuspidal homology
  vector<int> needed, freegens;
  long rk, denom1, denom2, dimension, denom3, ncusps;
  edge_relations ER;
  face_relations FR;
  int plusflag;
  Qideal N; // the level
  P1N P1;
  long ngens, nsymb, nap;
  vector<Quadprime> primelist; // list of primes in order, bad primes first

  ssubspace kern;
  smat tkernbas; // transpose of kernel(delta) basis
  vector<modsym> freemods;
  mat projcoord;
  long hmod; // if >0, failed to lift from modular linear algebra
             // so coord is modulo this

  homspace(const Qideal& I, int hp, int cuspid, int verb);

  void kernel_delta();          // computes ker(delta) for cuspidal homology
  void make_freemods();         // computes freemods and needed

  int check_conjugate(int verb=0);     // function to check consistency between this and conjugate level

  vec coords(int i) {return FR.coords(i);}
  int coords(const Quad& c, const Quad& d);
  int index(const Quad& c, const Quad& d) {return P1.index(c, d);}

   mat opmat(int, int dual=1, int verb=0);
   vec opmat_col(int, int j, int verb=0);
   mat opmat_cols(int, const vec& jlist, int verb=0);
   mat opmat_restricted(int i,const subspace& s, int dual, int verb=0);
  // versions returning an smat:
   smat s_opmat(int i,int dual, int verb=0);
   svec s_opmat_col(int i, int j, int verb=0);
   smat s_opmat_cols(int i, const vec& jlist, int verb=0);
   smat s_opmat_restricted(int i, const ssubspace& s, int dual,int verb=0);
   long matdim(void) {return dimension;}
   vector<long> eigrange(long i);

   long h1cuspdim() const {return dim(kern);}
   long h1dim() const {return dimension;}  // No confusion with subspace::dim
   long h1denom() const {return denom1;}
   long h1cdenom() const {return denom3;}
   long h1ncusps() const {return ncusps;}
   long h1hmod() const {return hmod;}

  vec chaincd(const Quad& c, const Quad& d, int type=0, int proj=0);
  vec chain(const Quad& a, const Quad& b, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one);
  vec chain(const RatQuad& r, int proj=0, const Quad& c=Quad::zero, const Quad& d=Quad::one)
  {return chain(r.num(), r.den(), proj, c, d);}
  vec chain(const RatQuad& alpha, const RatQuad& beta, int proj=0);
  vec chain(const modsym& m, int proj=0)
  {return chain(m.alpha(), m.beta(), proj);}

  vec applyop(const matop& mlist, const RatQuad& m, int proj=0);
  vec applyop(const matop& mlist, const modsym& m, int proj=0);

  mat calcop(const string opname, const matop& mlist, int dual=1, int display=0);
  vec calcop_col(const matop& mlist, int j);
  mat calcop_cols(const matop& mlist, const vec& jlist);
  smat s_calcop(const string opname, const matop& mlist, int dual, int display);
  smat s_calcop_cols(const matop& mlist, const vec& jlist);
  mat calcop_restricted(const string opname, const matop& mlist, const subspace& s, int dual, int display);
  smat s_calcop_restricted(const string opname, const matop& mlist, const ssubspace& s, int dual, int display);

public:
  string opname(const Quadprime& P);
  string opname(const Qideal& P);
   mat heckeop(Quadprime& P, int dual=1, int display=0);
   vec heckeop_col(Quadprime& P, int j, int display=0);
   mat heckeop_cols(Quadprime& P, const vec& jlist, int display=0);
   smat s_heckeop(Quadprime& P, int dual, int display);
   svec s_heckeop_col(Quadprime& P, int j, int display);
   smat s_heckeop_cols(Quadprime& P, const vec& jlist, int display);
   mat heckeop_restricted(Quadprime& P, const subspace& s, int dual, int display);
   smat s_heckeop_restricted(Quadprime& P, const ssubspace& s, int dual, int display);
   mat wop(Quadprime& Q, int dual=1, int display=0);
   mat fricke(int dual=1, int display=0);
//   mat conj(int display=0) const;
  vec maninvector(Quadprime& P, int proj=0);
  vec manintwist(const Quad& lambda, const vector<Quad>& res, vector<int> chitable, int proj=0);
  // no longer used, implementation commented out:
  // vec newhecke(const Quadprime& P, const Quad& n, const Quad& d);
};

// Each relation is a signed sum of edges (M)_alpha = {M(alpha},
// M(oo)} for M in the list mats and alpha=alphas[t] (when t>=0) or
// sigmas[-t] (when t<0), for t in the list types.  Here we check that such a
// relation holds identically in H_3 (not just modulo the congruence
// subgroup!)

// General case:
int check_rel(const vector<mat22>& mats, const vector<int>& types, const vector<int>& signs);
// Special case: all signs +1
int check_rel(const vector<mat22>& mats, const vector<int>& types);


vec reduce_modp(const vec& v, const scalar& p=DEFAULT_MODULUS);
mat reduce_modp(const mat& m, const scalar& p=DEFAULT_MODULUS);

// List of bad primes (dividing N) followed by good primes to length np:
vector<Quadprime> primelist(Qideal& N, int np);

#endif
