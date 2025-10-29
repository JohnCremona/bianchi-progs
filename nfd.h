// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "homspace.h"
#include "field.h"

class Newform;
class Newforms;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
private:
  Newforms* nf;    // pointer to "parent" class holding global info
  int index;       // index (starting from 1) of this newforms in the list of all
  int d;      // dim(S)
  Field* F;   // pointer to the (homological) Hecke field
  subspace S; // irreducible subspace of modular symbol space
  scalar denom_rel, denom_abs; // relative and absolute denominators of S
  vector<scalar> scales; // powers of denom_rel
  //  vector<scalar> contents;
  mat projcoord; // used to computed eigenvalues of any operator
  vector<int> epsvec;  // list of unramified quadratic character values (+1,-1) on S
  INT genus_char_disc; // associated discriminant factor (1 for trivial char)
  // book-keeping data for eigenvalue computations
  FieldModSq* Fmodsq;
  vector<long> genus_classes; // list of classes for which we have a nonzero eigenvalue
  vector<Qideal> genus_class_ideals; // list of squarefree ideals in these classes
  vector<Eigenvalue> genus_class_aP;  // list of eigenvalues of these ideals
  int genus_classes_filled;  // Set to 1 when all genus classes are
                             // filled, or when half are filled if we
                             // have detected self-twist
  // For each genus class c we count how many primes P in class c have
  // a(P)=0, to aid in detecting self-twist forms.
  vector<int> genus_class_trivial_counter;
  // list of possible self-twist discriminants, initially depends only
  // on the level but may be cut down later
  vector<INT> possible_self_twists;
  // map of eigenvalues for (good) primes, computed by geteigs()
  map<Quadprime, Eigenvalue> aPmap;
  // Fill dict aPmap of eigenvalues of first ntp good primes
  void compute_eigs(int ntp=10, int verbose=0);

public:
  // constructor from ambient Newforms using one irreducibel factor of char
  // poly of Newforms's T_mat
  Newform(Newforms* x, int ind, const ZZX& f, int verbose=0);
  // eigenvalue in F of a general principal operator on this:
  FieldElement eig(const matop& T);
  // eigenvalue of AutoHeckeOp(P) on this:
  FieldElement ap(Quadprime& P);
  // eigenvalue of a scalar operator
  ZZ eps(const matop& T);

  // eigenvalue of a (good) prime
  Eigenvalue eig(Quadprime& P);
  Field* field() const {return F;}
  string var() const {return F->var;}
  // output basis for the Homological Hecke field and character
  // If full, also output multiplicative basis for the full Hecke field
  void display(int full=0);
  int dimension() const {return d;}
  ZZX poly() const {return F->poly();}
  vector<int> character() const {return epsvec;}
  int trivial_char(); // 1 iff unramified quadratic character values (if any) are all +1
  ZZ basis_factor() const {return F->Bdet;}
  map<Quadprime, Eigenvalue> eigs(int ntp=10, int verbose=0)
  {
    compute_eigs(ntp, verbose);
    return aPmap;
  }
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
  friend class Field;
  friend class FieldElement;
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

  mat_m T_mat;  // matrix of splitting operator
  string T_name;  // name of splitting operator
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)
  // list of possible self-twist discriminants depending only on the level
  vector<INT> possible_self_twists;

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

  mat_m heckeop(Quadprime& P, int cuspidal=0, int dual=0); // not const as may add info into N
  mat_m heckeop(const matop& T, int cuspidal=0, int dual=0) const;
  mat_m heckeop(const gmatop& T, int cuspidal=0, int dual=0) const;

  int ok() const {return split_ok;}
  int nforms() const {return newforms.size();}
  string splitopname() const {return T_name;}
  vector<int> dimensions() const;
  // output basis for the Hecke field and character of all newforms
  void display_newforms(int triv_char_only=0, int full=0) const;
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}
};

#endif
