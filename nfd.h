// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar class (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "homspace.h"

class Newforms;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
private:
  Newforms* nf;    // pointer to "parent" class holding global info
  int d;      // dim(S)
  ZZX minpoly;  // irredicible poly defining S of degree d
  subspace S; // irreducible subspace of modular symbol space
  mat A;        // matrix of T restricted to S, with min.poly. minpoly
  scalar denom_rel, denom_abs; // relative and absolute denominators of S
  vector<scalar> scales; // powers of denom_rel
  vector<scalar> contents;
  mat projcoord; // used to computed eigenvalues of any operator
  vector<scalar> epsvec; // list of unramified quadratic character values on S
  scalar nfbasis_num, nfbasis_den; // num&denom of basis scale factor
  mat nfbasis;   // columns give newform basis in terms of powers basis
public:
  // constructor from ambient Newforms using one irreducibel factor of char
  // poly of Newforms's T_mat
  Newform(Newforms* x, const ZZX& f, int verbose=0);
  // eigenvalue (as coords w.r.t. basis) of a general operator on this:
  vec eig(const matop& T) const;
  // eigenvalue of AutoHeckeOp(P) on this:
  vec ap(Quadprime& P) const;
  // eigenvalue of a scalar operator
  scalar eps(const matop& T) const;

  // output basis for the Hecke field and character
  void display_basis(int j) const;
  int dimension() const {return d;}
  ZZX poly() const {return minpoly;}
  vector<scalar> character() const {return epsvec;}
};

// class for the collection of all d-dimensional newforms
class Newforms {
  friend class Newform;
private:
  int verbose;

  Qideal N; // the level
  vector<Qideal> Ndivs; // divisors of N
  vector<matop> eps_ops; // list of unram quad chars

  homspace* H1;  // the ambient modular symbol space at level N
  int cdimH, dimH;
  scalar hmod, dH;
  vector<scalar> Hscales;

  mat T_mat;  // matrix of splitting operator
  string T_name;  // name of splitting operator
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)
  vector<Newform> newforms; // the newforms

  // Internal methods, called by constructor

  //void find_T_manual(); // compute T (via prompts)
  //void factor_T();

  // If simple find T_P whose char poly on the newspace is squarefree
  // and coprime to its char poly on the oldspace, trying all good P
  // with N(P)<=maxnormP.  If not simple, try similar where T is a
  // linear combination of T_P.  Set split_ok=1 if successful else 0.
  void find_T(int simple=1, INT maxnormP=INT(0));

public:
  Newforms(void) {;}
  // constructor from a homspace, looking for a splitting operator
  // from primes up to given norm
  Newforms(homspace* h1, const INT& maxnormP, int verb=1);
  int split_ok; // records whether the constructor was able to find a splitting operator
  mat heckeop(Quadprime& P, int cuspidal=0, int dual=0);
  mat heckeop(const matop& T, int cuspidal=0, int dual=0) const;
  mat heckeop(const gmatop& T, int cuspidal=0, int dual=0) const;
  vector<vec> ap(Quadprime& P);
  vector<vec> eig(const matop& T) const;
  //  vector<scalar> eps(const matop& T) const; // T should be a scalar operator

  int ok() const {return split_ok;}
  int nforms() const {return newforms.size();}
  string splitopname() const {return T_name;}
  vector<int> dimensions() const;
  // output basis for the Hecke field and character of all newforms
  void display_bases() const;
};


#endif
