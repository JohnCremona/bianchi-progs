// File NFD.H: class for newforms of any dimension
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "homspace.h"
#include "field.h"

class Newform;
class Newspace;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
  friend class Newspace;
private:
  Newspace* nf;    // pointer to "parent" class holding global info
  int index;       // index (starting from 1) of this newforms in the list of all
  string lab;    //
  int d;      // dim(S)
  Field* F;   // pointer to the (homological) Hecke field
  subspace S; // irreducible subspace of modular symbol space
  scalar denom_rel, denom_abs; // relative and absolute denominators of S
  vector<scalar> scales; // powers of denom_rel
  //  vector<scalar> contents;
  mat projcoord; // used to computed eigenvalues of any operator
  vector<int> epsvec;  // list of unramified quadratic character values (+1,-1) on S
  int triv_char;  // 1 iff all epsvec values are +1, else 0
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

  int self_twist_flag;
  INT CMD;            // =D if this is self-twist by unramified disc D dividing Quad::disc, else 0
  int bc;  // base-change code (0 for no, 1 for yes, -1 for unknown)
  int bct; // base-change twist code (0 for no, 1 for yes including 1 when bc==1, -1 for unknown)
  int cm; // CM code (see defn in newforms.h, but only 1 (not set) for now)

  // Dict of eigenvalues of principal operators (the key includes an
  // int for sorting, othewise they get sorted in alphabetical order
  // of opname):
  map<pair<int,string>, Eigenvalue> eigmap;
  // Dict of T(P) eigenvalues of good primes P:
  map<Quadprime, Eigenvalue> aPmap;
  // Dict of W(Q) eigenvalues of bad primes Q:
  map<Quadprime, Eigenvalue> eQmap;
  // Fill dict aPmap of eigenvalues of first ntp good primes
  void compute_eigs(int ntp=10, int verbose=0);
  // Fill dict eigmap of eigenvalues of first ntp principal operators
  Eigenvalue compute_one_principal_eig(int, const matop& T, int store=0, int verbose=0);
  void compute_principal_eigs(int ntp=10, int verbose=0);
  // Fill dict aPmap of eigenvalues of first ntp good primes, trivial char case only
  void compute_eigs_triv_char(int ntp=10, int verbose=0);
  // Fill dict aPmap of eigenvalues of first ntp good primes, class group C4 only
  void compute_eigs_C4(int ntp=10, int verbose=0);
  // Fill dict eQmap *after* aPmap, if triv_char
  void compute_AL_eigs(int verbose=0);
  // Assuming aPmap filled, set the bc and bct flags
  void check_base_change(void);

public:
  // constructor from ambient Newspace using one irreducibel factor of char
  // poly of Newspace's T_mat
  Newform(Newspace* x, int ind, const ZZX& f, int verbose=0);

  int get_index() const { return index;}
  void set_index(int i) {index = i; lab = codeletter(i-1); F->set_var(lab);}

  // Functions for computing eigenvalues of principal operators:

  // eigenvalue of a general principal operator:
  FieldElement eig(const matop& T);
  // eigenvalue of AutoHeckeOp(P):
  FieldElement ap(Quadprime& P);
  // eigenvalue of a scalar operator
  ZZ eps(const matop& T);

  // eigenvalue of a (good) prime from aPmap if P is in there;
  // otherwise either raise an error (if stored_only=1) or (not yet
  // implemented) compute it.
  Eigenvalue eig(Quadprime& P, int stored_only=1);

  // Principal eigenvalue of a (good) prime P if P has trivial genus
  // class, or P^2 otherwise, from aPmap.  Assuming trivial character
  // this will be the eigenvalue of AutHeckeOp(P):
  FieldElement eigPorP2(Quadprime& P);
  // Principal eigenvalue of a linear combination of the above:
  FieldElement eig_lin_comb(vector<Quadprime>& Plist, vector<scalar>& coeffs);
  // Characteristic polynomial of such a linear combination:
  ZZX char_pol_lin_comb(vector<Quadprime>& Plist, vector<scalar>& coeffs);

  Field* field() const {return F;}
  string label() const {return lab;}
  // output basis for the Homological Hecke field and character
  // If full, also output multiplicative basis for the full Hecke field
  void display(int full=0);
  int dimension(int full=1) const
  {
    if (full)
      return d<<Fmodsq->rank();
    else
      return d;
  }
  ZZX poly() const {return F->poly();}
  vector<int> character() const {return epsvec;}
  int is_char_trivial() const {return triv_char;}
  int is_self_twist() const {return self_twist_flag;}
  INT self_twist_discriminant() const {return CMD;}
  // return +1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know
  int is_base_change(void)
  {
    // try to determine bc, bct if not yet done:
    if (bc==-1) check_base_change();
    if (bc==1) return 1;
    else if (bct==1) return -1;
    else if (bct==0) return 0;
    else return 2;
  }
  int cm_code() const {return cm;}
  ZZ basis_factor() const {return F->Bdet;}
  map<Quadprime, Eigenvalue> aPeigs(int ntp=10, int verbose=0)
  {
    if (aPmap.empty())
      compute_eigs(ntp, verbose);
    return aPmap;
  }
  map<Quadprime, Eigenvalue> ALeigs(int verbose=0)
  {
    if (eQmap.empty() && triv_char)
      compute_AL_eigs(verbose);
    return eQmap;
  }
  map<pair<int,string>, Eigenvalue> principal_eigs(int nap=10, int verbose=0)
  {
    compute_principal_eigs(nap, verbose);
    return eigmap;
  }
  // filename for this Newform
  string filename() const;
  void output_to_file() const;
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
class Newspace {
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
  Newspace(void) {;}
  // constructor from a homspace, looking for a splitting operator
  // using linear combinations of up to maxnp primes, coefficients up
  // to maxc
  Newspace(homspace* h1, int maxnp, int maxc, int verb=1);
  int split_ok; // records whether the constructor was able to find a splitting operator

  mat_m heckeop(Quadprime& P, int cuspidal=0, int dual=0); // not const as may add info into N
  mat_m heckeop(const matop& T, int cuspidal=0, int dual=0) const;
  mat_m heckeop(const gmatop& T, int cuspidal=0, int dual=0) const;

  int ok() const {return split_ok;}
  int nforms() const {return newforms.size();}
  string splitopname() const {return T_name;}
  vector<int> dimensions(int full=0) const;
  // output basis for the Hecke field and character of all newforms
  void display_newforms(int triv_char_only=0, int full=0) const;
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}
  // filename for Newspace
  string filename();
  // output data for this Newspace and each Newform
  void output_to_file();
};


#endif
