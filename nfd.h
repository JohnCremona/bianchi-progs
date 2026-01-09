// File NFD.H: class Newform for newforms of any dimension
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

  int self_twist_flag; // -1 for unknown, 0 for no, 1 for yes
  INT CMD;            // =D if this is self-twist by unramified disc D<0 dividing Quad::disc, else 0
  int bc;  // base-change code (0 for no, 1 for yes, -1 for unknown)
  int bct; // base-change twist code (0 for no, 1 for yes including 1 when bc==1, -1 for unknown)
  int cm; // CM code (see defn in newforms.h, but only 1 (not set) for now)
  int sfe; // Sign of functional equation (minus product of AL eigs), 0 if not known

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
  // constructor from ambient Newspace using one irreducible factor of char
  // poly of Newspace's T_mat
  Newform(Newspace* x, int ind, const ZZX& f, int verbose=0);
  // constructor from ambient Newspace (read from file)
  Newform(Newspace* x, int i, int verbose=0);

  int get_index() const { return index;}
  void set_index(int i) {index = i; lab = codeletter(i-1); F->set_var(lab);}

  // Functions for computing eigenvalues of principal operators:

  // eigenvalue of a general principal operator:
  FieldElement eig(const matop& T);
  // eigenvalue of AutoHeckeOp(P):
  FieldElement ap(const Quadprime& P);
  // eigenvalue of a scalar operator
  ZZ eps(const matop& T);

  // eigenvalue of a (good) prime from aPmap if P is in there;
  // otherwise either raise an error (if stored_only=1) or (not yet
  // implemented) compute it.
  Eigenvalue eig(const Quadprime& P, int stored_only=1);

  // Principal eigenvalue of AutoHeckeOp(P) for a good prime P, from
  // stored aP in aPmap.  Only implemented for trivial character
  // (where this is the eigenvalue of P or P^2) or C4 class group.
  FieldElement eigPauto(const Quadprime& P, int verb=0);
  // Principal eigenvalue of a linear combination of the above:
  FieldElement eig_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs, int verb=0);
  // Characteristic polynomial of such a linear combination:
  ZZX char_pol_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs, int verb=0);

  Field* field() const {return F;}
  string label_suffix() const {return lab;}
  string short_label() const; // level_label-suffix
  string long_label() const;  // field_label-level_label-suffix
  string conj_label() const; // conj-level_label-suffix
  string long_conj_label() const;  // field_label-conj-level_label-suffix
  // Output basis for the Homological Hecke field and character
  // If class number even, also output multiplicative basis for the full Hecke field
  // Optionally aP and AL data too
  void display(int aP=0, int AL=0, int principal_eigs=0) const;
  // Display aP data (trivial char or C4 fields)
  void display_aP() const;
  // Display A-L eigenvalues (trivial char or C4 fields)
  void display_AL() const;
  // Display principal eigenvalues
  void display_principal_eigs() const;

  int dimension(int full=1) const;
  ZZX poly() const {return F->poly();}
  vector<int> character() const {return epsvec;}
  int is_char_trivial() const {return triv_char;}
  int is_self_twist() const {return self_twist_flag;}
  INT self_twist_discriminant() const {return CMD;}
  // return base-change code (+1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know)
  int base_change_code(void) const;
  int cm_code() const {return cm;}
  ZZ basis_factor() const {return F->Bdet;}
  // Compute aPmap if empty and return it
  map<Quadprime, Eigenvalue> aPeigs(int ntp, int verbose=0);
  // Compute eQmap if empty and return it
  map<Quadprime, Eigenvalue> ALeigs(int verbose=0);
  // Compute eigmap (principal eigs) if empty and return it
  map<pair<int,string>, Eigenvalue> principal_eigs(int nap, int verbose=0);

  // NB We only implement file output for newforms with trivial character
  // Filename for this Newform (or conjugate):
  string filename(int conj=0) const;
  // Output newform data (or data for the conjugate newform):
  void output_to_file(int conj=0) const;
  // Input newform data. Returns 0 if data not available, else 1.
  int input_from_file(int verb=0);
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
  string lab; // the level's label
  vector<Qideal> Ndivs; // divisors of N
  vector<Quadprime> badprimes; // prime divisrs of N
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

  map<string, Newspace> oldspaces; // keys are labels of all proper divisors
  int make_oldspaces(); // reads Newspaces for proper divisors from file, 1 iff success
  // return full and new char polys for a linear combo of ops using old Newspace data from files
  pair<ZZX,ZZX> full_and_new_polys(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                                   const gmatop &T);
  // Return true iff this combo of ops has squarefree new poly coprime to its old poly
  int valid_splitting_combo(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                            const gmatop &T, ZZX& f_new);

  // Find a linear combination of up to maxnp operators (T_{A,A} or
  // T_P) with coefficients up to maxc, whose char poly on the
  // newspace is squarefree and coprime to its char poly on the
  // oldspace.  Set split_ok=1 if successful else 0.
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
  Qideal level() const {return N;}
  string short_label() const {return lab;}
  string long_label() const {return field_label() + "-" + lab;}
  string conj_label() const {Qideal Nbar = N.conj();return ideal_label(Nbar);}
  string long_conj_label() const {return field_label() + "-" + conj_label();}

  string splitopname() const {return T_name;}
  vector<int> dimensions(int full=0) const;

  // output all newforms: Dimension, Character, Hecke field; optionally aP and AL data
  void display_newforms(int aP=0, int AL=0, int principal_eigs=0, int triv_char_only=0) const;
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}

  // filename for Newspace (or conjugate)
  string filename(int conj=0) const;
  // output data for this Newspace (or conjugate) and each Newform
  void output_to_file(int conj=0) const;
  // Input Newspace data and newform data for each newform. Returns 0 if data missing, else 1.
  int input_from_file(const Qideal& level, int verb=0);
};

// test whether field's class group is C4
inline int is_C4()
{
  return (Quad::class_number==4) && (Quad::class_group_2_rank==1); // C4
}

// For class group C4 only (so far)
// return v where v[i] is the index of ideal class c^i for one generator class c
vector<int> C4classes();

#endif
