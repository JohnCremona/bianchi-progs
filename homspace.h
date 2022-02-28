// FILE HOMSPACE.H: Declaration of class homspace

#if     !defined(_HOMSPACE_H)
#define _HOMSPACE_H      1       //flags that this file has been included

#include "cusp.h"
#include "face_relations.h"
#include "hecke.h"

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
  long ngens, nsymb, nap, nwq;

  ssubspace kern;
  smat tkernbas; // transpose of kernel(delta) basis
  vector<modsym> freemods;
  mat projcoord;
  long characteristic; // =0 or prime p
  long hmod; // if >0, did not lift from modular linear algebra
             // so coord is modulo this

  homspace(const Qideal& I, int hp, int cuspid, int verb=0, long ch=0);

  void kernel_delta();          // computes ker(delta) for cuspidal homology
  void make_freemods();         // computes freemods and needed

  int check_conjugate(int verb=0);     // function to check consistency between this and conjugate level

  vec coords(int i) {return FR.coords(i);}
  int coords(const Quad& c, const Quad& d);
  int index(const Quad& c, const Quad& d) {return P1.index(c, d);}

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

  vec applyop(const matop& T, const RatQuad& m, int proj=0);
  vec applyop(const matop& T, const modsym& m, int proj=0);

  mat calcop(const matop& T, int dual=1, int display=0);
  vec calcop_col(const matop& T, int j, int verb=0)  {return applyop(T,freemods[j-1]);}
  mat calcop_cols(const matop& T, const vec& jlist, int verb=0);
  mat calcop_restricted(const matop& T, const subspace& s, int dual, int display);
  smat s_calcop(const matop& T, int dual, int display);
  smat s_calcop_cols(const matop& T, const vec& jlist, int verb=0);
  smat s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display);

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

#endif
