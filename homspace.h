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
  long ngens, nsymb, nap, n2r, nwq;
  vector<Quadprime> badprimes; // list of bad primes Q with square ideal class or even exponent
  vector<Quadprime> goodprimes;  // good primes in order
  vector<Qideal> nulist; // list of ideals coprime to level generating 2-torsion in class group

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

  // For the automatic finding of 1-dimensional egenspaces we need to
  // interface with eclib's splitter class, which wants to know thw
  // i'th operator matix and its possible eigenvalues, for i=0,1,2,...

  // the list of matrices defining the i'th operator:
  matop h1matop(int);
  // the list of possible (integer) eigenvalues for the i'th operator:
  vector<long> eigrange(long i);

  // the next 9 functions are required by the splitter_base class (via
  // functions of the same name in the newforms class which is derived
  // from the splitter_base class)

  // mat/vec versions:
  mat opmat(int, int dual=1, int verb=0);
  vec opmat_col(int, int j, int verb=0);
  mat opmat_cols(int, const vec& jlist, int verb=0);
  mat opmat_restricted(int i,const subspace& s, int dual, int verb=0);
  // smat/svec versions:
  smat s_opmat(int i,int dual, int verb=0);
  svec s_opmat_col(int i, int j, int verb=0);
  smat s_opmat_cols(int i, const vec& jlist, int verb=0);
  smat s_opmat_restricted(int i, const ssubspace& s, int dual,int verb=0);
  long matdim(void) {return dimension;}

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
  vec calcop_col(const matop& T, int j)  {return applyop(T,freemods[j-1]);}
  mat calcop_cols(const matop& T, const vec& jlist);
  mat calcop_restricted(const matop& T, const subspace& s, int dual, int display);
  smat s_calcop(const matop& T, int dual, int display);
  smat s_calcop_cols(const matop& T, const vec& jlist);
  smat s_calcop_restricted(const matop& T, const ssubspace& s, int dual, int display);

public:
  mat heckeop(Quadprime& P, int dual=1, int display=0)
  {
    return calcop(HeckeOp(P,N), dual, display);
  }

  vec heckeop_col(Quadprime& P, int j, int display=0)
  {
    return calcop_col(HeckeOp(P,N), j);
  }

  mat heckeop_cols(Quadprime& P, const vec& jlist, int display=0)
  {
    return calcop_cols(HeckeOp(P,N), jlist);
  }

  smat s_heckeop(Quadprime& P, int dual, int display)
  {
    return s_calcop(HeckeOp(P,N), dual, display);
  }

  svec s_heckeop_col(Quadprime& P, int j, int display)
  {
    return calcop_col(HeckeOp(P,N), j);
  }

  smat s_heckeop_cols(Quadprime& P, const vec& jlist, int display)
  {
    return s_calcop_cols(HeckeOp(P,N), jlist);
  }

  mat heckeop_restricted(Quadprime& P, const subspace& s, int dual, int display)
  {
    return calcop_restricted(HeckeOp(P,N), s, dual, display);
  }

  smat s_heckeop_restricted(Quadprime& P, const ssubspace& s, int dual, int display)
  {
    return s_calcop_restricted(HeckeOp(P,N), s, dual, display);
  }

  mat wop(Quadprime& Q, int dual=1, int display=0)
  {
    return calcop(AtkinLehnerOp(Q,N), dual,display);
  }

  mat fricke(int dual=1, int display=0)
  {
    return calcop(FrickeOp(N), dual,display);
  }

  // unramified character (A coprime to the level and A^2 principal)
  mat nu_op(Qideal& A, int dual=1, int display=0)
  {
    return calcop(CharOp(A, N), dual,display);
  }
  // as previous, for the i'th such involution (for 0<=i<Quad::class_group_2_rank)
  mat nu(int i, int dual=1, int display=0)
  {
    return nu_op(nulist[i], dual, display);
  }

  // T_{P^2} when P^2 principal
  mat hecke_sq_op(Quadprime& P, int dual=1, int display=0)
  {
    return calcop(HeckeSqOp(P, N), dual,display);
  }

  // (T_P)^2 when P^2 principal, using  (T_P)^2 = T_{P^2} + N(P)*T_{P,P}
  mat hecke_op_sq(Quadprime& P, int dual=1, int display=0)
  {
    return calcop(HeckeSqOp(P, N), dual,display)+ I2long(P.norm())*nu_op(P, dual, display);
  }

  // T_{PQ} for P!=Q when PQ principal
  mat hecke_pq_op(Quadprime& P, Quadprime& Q, int dual=1, int display=0)
  {
    return calcop(HeckePQOp(P, Q, N), dual,display);
  }

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

// List of bad primes (dividing N) followed by good primes to length
// at least np, making sure that the list includes at least one good
// principal prime.  iP0 is set to the index in the list of the first
// good principal prime.
// If p (default 0) is nonzero, omit bad primes and primes dividing P
vector<Quadprime> make_primelist(Qideal& N, int np, int& iP0, int p=0);

#endif
