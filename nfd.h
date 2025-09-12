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
  homspace* H1;  // the ambient modular symbol space at level N
  mat_m T;  // matrix of an operator on H1

  mat_m W,Winv,Winv_scaled;
  int cdimH, dimH;
  scalar coord_fac;
  bigint hmod, Wdetnum, Wdetdenom;
  vector<bigint> minpol;  // min poly of alpha, field generator
  int verbose;
public:
  ZZX f;    // a multiplicity-1 irreducible factor of charpoly(T)
  subspace_m S;  // the irreducible subspace, ker(f(T))
  mat_m A;  // matrix of T restricted to S
  vector<bigint> Hscales, Sscales;
  bigint dH, dS, dHS;
  int dimS;
  nfd(void) {;}
  nfd(homspace* h1, int verb=1);
  void display(void) const;
  mat_m heckeop(Quadprime& P);
  mat_m heckeop_S(Quadprime& P);
  vec_m ap(Quadprime& P);
  void make_T(); // compute T (via prompts)
  void make_S(); // compute S, A, f (via prompts)
};

extern vector<long> class_number_one_fields; // defined in quads.cc
int is_class_number_one(long d);

extern map<Qideal,homspace*> H1_dict;
homspace* get_homspace(const Qideal& N, scalar mod);

extern map<pair<Qideal,Quadprime>, ZZX> full_poly_dict;
ZZX get_full_poly(const Qideal& N,  Quadprime& P, const scalar& mod);

extern map<pair<Qideal,Quadprime>, ZZX> new_poly_dict;
ZZX get_new_poly(Qideal& N, Quadprime& P, const scalar& mod);

#endif
