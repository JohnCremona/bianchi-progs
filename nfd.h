// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar class (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "homspace.h"

class nfd {
private:
  Qideal N; // the level
  vector<Qideal> Ndivs; // divisors of N
  homspace* H1;  // the ambient modular symbol space at level N

  int cdimH, dimH;
  scalar hmod, dH;
  int verbose;
  vector<scalar> Hscales;

  mat T;  // matrix of an operator on H1

  // data per irreducible component
  vector<subspace> S;  // the irreducible subspaces, ker(f(T)) for f in factors
  vector<scalar> dS, dHS; // relative and absolute denominator of S
  vector<vector<scalar>> Sscales, Scontents;
  vector<mat> A;  // matrices of T restricted to S
  vector<mat> W,Winv,Winv_scaled;
  vector<scalar> Wdetnum, Wdetdenom;
  vector<mat> projcoord;
public:
  nfd(void) {;}
  nfd(homspace* h1, int verb=1);
  void display(void) const;
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)
  int nfactors;        // the number of them
  vector<int>dimS;     // their degrees
  void display_basis(int j) const; // output basis info for subspace j (1<=j<=nfactors)

  mat heckeop(Quadprime& P, int cuspidal=0, int dual=0);
  vector<vec> ap(Quadprime& P);
  void find_T(); // compute T (via prompts)
  // compute T=T_P, trying all good P with N(P)<=maxnormP
  int find_T_auto(INT maxnormP, Quadprime& P0, int verb);
  void factor_T();
  void make_irreducible_subspaces(); // compute S, A, f for all factors
};

extern vector<long> class_number_one_fields; // defined in quads.cc
int is_class_number_one(long d);

#endif
