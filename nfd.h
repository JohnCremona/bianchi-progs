// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////
//
// Adapted from the similar class (over Q) in eclib
//
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "matprocs.h"
#include "homspace.h"

class HeckeFieldElement;

class HeckeField {
  friend class HeckeFieldElement;
private:
  int d;        // degree
  ZZX minpoly;  // irredicible poly of degree d
  scalar denom; // minpoly is the (integral) min poly of A/denom
  mat A;        // dxd matrix with scaled min.poly. minpoly
  mat C;        // dxd companion matrix with min.poly. minpoly
  mat B, Binv;  // Binv*B = Bdet*I
  scalar Bdet;  // Binv*A*B = Bdet*denom*C
  scalar Bfactor; // basis scale factor: cols of Binv are Bfactor * coeffs of basis w.r.t. a-powers
public:
  //HeckeField(const ZZX& p);
  HeckeField(); // defaults to Q
  HeckeField(const mat& m, const scalar& denom = scalar(1), int verb=0);
  int degree() const {return d;}
  ZZX poly() const {return minpoly;}
  mat basis() const {return Binv;} // columns are Bfactor * coeffs of basis w.r.t. a-powers
  scalar basis_factor() const {return Bfactor;}
  mat inv_basis() const {return B;} // columns are coeffs of a-powers w.r.t. basis
  void display(ostream&s = cout) const;
};

class Newforms;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
private:
  Newforms* nf;    // pointer to "parent" class holding global info
  int d;      // dim(S)
  HeckeField F;
  subspace S; // irreducible subspace of modular symbol space
  scalar denom_rel, denom_abs; // relative and absolute denominators of S
  vector<scalar> scales; // powers of denom_rel
  //  vector<scalar> contents;
  mat projcoord; // used to computed eigenvalues of any operator
  vector<int> epsvec;  // list of unramified quadratic character values (+1,-1) on S
  INT genus_char_disc; // associated discriminant factor (1 for trivial char)
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
  void display(int j) const; // j is the index in the list of all newforms
  int dimension() const {return d;}
  ZZX poly() const {return F.poly();}
  vector<int> character() const {return epsvec;}
  int trivial_char(); // 1 iff unramified quadratic character values (if any) are all +1
};

// function to sort newforms of the same level, by (1) character
// values (reverse lexicographically so trivial char is first), (2)
// dimension, (3) min poly (so dimension 1 forms with the same
// character are ordered only by the order of the degree 1 factors of
// the splitting polynomial).

struct newform_comparison {
  bool operator()(Newform& f1, Newform& f2)
  {
    // first sort by character
    vector<int> char1 = f1.character(), char2 = f2.character();
    bool t = char1>char2;  // true if f1 has 'earlier' char
    if(t) return 1;
    t = char1<char2; // true if f1 has 'later' char
    if(t) return 0;

    // Now sharacters are the same,  sort by dimension
    int s = f1.dimension() - f2.dimension();
    if(s) return (s<0); // true if f1 has smaller dimension

    // then sort by min poly (smaller degree comes before larger)
    ZZX pol1 = f1.poly(), pol2 = f2.poly();
    return poly_cmp(pol1, pol2);
  }
};

extern newform_comparison newform_cmp;

// class for the collection of all d-dimensional newforms
class Newforms {
  friend class Newform;
  friend class HeckeField;
  friend class HeckeFieldElement;
private:
  int verbose;

  Qideal N; // the level
  vector<Qideal> Ndivs; // divisors of N
  vector<Qideal> t2ideals; // list of ideals coprime to level generating 2-torsion in class group
  vector<matop> eps_ops; // list of T(A,A) operators for A in t2ideals

  homspace* H1;  // the ambient modular symbol space at level N
  int cdimH, dimH;
  scalar hmod, dH;
  vector<scalar> Hscales;

  mat T_mat;  // matrix of splitting operator
  string T_name;  // name of splitting operator
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)

  // Internal methods, called by constructor

  //void find_T_manual(); // compute T (via prompts)
  //void factor_T();

  // If maxc=0, find T_P whose char poly on the newspace is squarefree
  // and coprime to its char poly on the oldspace, trying the first
  // maxnp good primes.  If maxc>0, try similar where T is a linear
  // combination of up to maxnp T_P with coefficients up to maxc.  Set
  // split_ok=1 if successful else 0.
  void find_T(int maxnp, int maxc);

public:
  vector<Newform> newforms; // the newforms
  Newforms(void) {;}
  // constructor from a homspace, looking for a splitting operator
  // using linear combinations of up to maxnp primes, coefficients up
  // to maxc
  Newforms(homspace* h1, int maxnp, int maxc, int verb=1);
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
  void display_newforms(int triv_char_only=0) const;
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}
};

// same as m.output(cout) except no newlines between rows
void output_flat_matrix(const mat& m, ostream&s = cout);

#endif
