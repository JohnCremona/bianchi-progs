// File NFD.H: class Newform for newforms of any dimension
//////////////////////////////////////////////////////////////////////////

#ifndef _BIANCHI_NFD_H
#define _BIANCHI_NFD_H      1

#include "homspace.h"
#include "field.h"

// In principle we can compute newforms with arbitrary unramified
// character, via homological newforms with arbitrary quadratic
// character on Cl[2]. At present the implementation is only complete
// for newforms with trivial character (so, all newforms over fields
// of odd class number) and newforms with nontrivial character over
// fields with class group C4.  With the current implementation of the
// Eigenvalue class this could be extended to all class-groups whose
// 2-primary component has exponent 4, but for complete generality the
// Eigenvalue class would need to be extended to include higher order
// 2-power roots of unity.

// test whether field's class group is C4
inline int is_C4()
{
  return (Quad::class_number==4) && (Quad::class_group_2_rank==1); // C4
}

class Newform;
class Newspace;

// class for a d-dimensional newform, defined by an irreducible factor
// of the characteristic polynomial of some splitting operator T
class Newform {
  friend class Newspace;
private:
  Newspace* nf;    // pointer to "parent" class holding global info
  int index;       // index (starting from 1) of this newforms in the list of all
  string lab;      // label suffix (a,b,c,...)
  int d;      // dim(S)
  Field* F0;   // pointer to the (homological) Hecke field (original)
  Field* F;    // pointer to the (homological) Hecke field (reduced)
  FieldIso Fiso; // isomorphism from F0 to F (possibly identity)
  FieldModSq* Fmodsq; // relative full Hecke field as extension of F
  Field* Fabs;   // pointer to absolute full Hecke field
  FieldIso abs_emb; // isomorphism from F to Fabs (possibly identity)
  vector<FieldElement> im_gens;
  subspace S; // irreducible subspace of modular symbol space
  scalar denom_abs; // absolute denominator of S
  mat projcoord; // used to computed eigenvalues of any operator
  vector<int> epsvec;  // list of unramified quadratic character values (+1,-1) on S
  int triv_char;  // 1 iff all epsvec values are +1, else 0
  INT genus_char_disc; // associated discriminant factor (1 for trivial char)

  // book-keeping data for eigenvalue computations
  vector<long> genus_classes_no_ext; // list of classes for which we have a nonzero eigenvalue in F itself
  vector<Qideal> genus_class_no_ext_ideals; // list of ideals in these classes
  vector<long> genus_classes_nonzero; // list of classes for which we have a nonzero eigenvalue
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
  // max norm(P) for P in aPmap:
  INT maxP;
  // Dict of W(Q) eigenvalues in {+1,-1} of bad primes Q, triv char only:
  map<Quadprime, int> eQmap;
  // Dict of coefficients in Fabs of integral ideals M. Trivial
  // character only.
  map<Qideal, FieldElement> aMmap;
  vector<ZZ> trace_list; // list of traces of sorted integral ideals

  // Fill dict eigmap of eigenvalues of first ntp principal operators
  Eigenvalue compute_one_principal_eig(int, const matop& T, int store=0, int verbose=0);
  void compute_principal_eigs(int ntp=10, int verbose=0);
  // Fill dict aPmap of eigenvalues of first ntp good primes, trivial char case only
  void compute_eigs_triv_char(int ntp=10, int verbose=0);
  // Fill dict aPmap of eigenvalues of first ntp good primes, class group C4 only
  void compute_eigs_C4(int ntp=10, int verbose=0);
  // Fill dict eQmap, if triv_char; if aPmap not already filled, first
  // compute ntp aP
  void compute_AL_eigs(int ntp=10, int verbose=0);

public:
  // Fill dict aPmap of eigenvalues of first ntp good primes; put max norm(P) into maxP
  void compute_eigs(int ntp=10, int verbose=0);
  // Assuming aPmap filled, fill aMmap (Fourier coefficients), and
  // trace_list (ordered list of traces) using all P in aPmap if aPmap
  // already filled, else first compute ntp aP.
  void compute_coefficients(int ntp=10, int verbose=0);
  // Assuming aPmap filled, set the bc and bct flags
  void check_base_change(void);

  // constructor from ambient Newspace using one irreducible factor of char
  // poly of Newspace's T_mat
  Newform(Newspace* x, int ind, const ZZX& f, int verbose=0);
  // constructor from ambient Newspace (read from file)
  Newform(Newspace* x, int i, int verbose=0);

  int get_index() const { return index;}
  void set_index(int i) {index = i; lab = codeletter(i-1); F0->set_var(lab+string("0")); F->set_var(lab);}

  // Functions for computing eigenvalues of principal operators:

  // eigenvalue of a general principal operator:
  FieldElement eig(const matop& T);
  // eigenvalue of AutoHeckeOp(P):
  FieldElement ap(const Quadprime& P);
  // eigenvalue +-1 of a scalar involution operator
  int eps(const matop& T);

  // eigenvalue of a (good) prime from aPmap if P is in there;
  // otherwise either raise an error (if stored_only=1) or (not yet
  // implemented) compute it.
  Eigenvalue eig(const Quadprime& P, int stored_only=1);

  // coefficient in Fabs of integral ideal M from aMmap or computed
  // (and stored in aMmap) using multiplicative relations. Trivial
  // character only.
  FieldElement aM(Qideal& M); // not const Qideal& as we factor it

  // Principal eigenvalue of AutoHeckeOp(P) for a good prime P, from
  // stored aP in aPmap.  Only implemented for trivial character
  // (where this is the eigenvalue of P or P^2) or C4 class group.
  // 'biglevel' is a multiple of the current level, auxiliary ideals A
  // must be coprime to this, not just to the current level
  FieldElement eigPauto(Quadprime& P, const Qideal& biglevel, int verb=0);
  // Principal eigenvalue of a linear combination of the above:
  FieldElement eig_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                            const Qideal& biglevel, int verb=0);
  // Characteristic polynomial of such a linear combination:
  ZZX char_pol_lin_comb(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                        const Qideal& biglevel, int verb=0);

  Field* field(int original=0) const {return (original? F0: F);}
  // Return the degree of the principal or full Hecke field
  int dimension(int full=1) const {return (full? d<<Fmodsq->rank() : d);}
  ZZX poly(int original=0) const {return (original? F0->poly(): F->poly());}
  string label_suffix() const {return lab;}
  string short_label() const; // level_label-suffix
  string long_label() const;  // field_label-level_label-suffix
  string conj_label() const; // conj-level_label-suffix
  string long_conj_label() const;  // field_label-conj-level_label-suffix

  // Output basis for the Homological Hecke field and character
  // If class number even, also output multiplicative basis for the full Hecke field
  // Optionally aP and AL data too
  void display(int aP=0, int AL=0, int principal_eigs=0, int traces=0) const;
  // Display aP data (trivial char or C4 fields)
  void display_aP() const;
  // Display AL eigenvalues (trivial char or C4 fields)
  void display_AL() const;
  // Display principal eigenvalues
  void display_principal_eigs() const;

  vector<int> character() const {return epsvec;}
  int is_char_trivial() const {return triv_char;}
  int is_self_twist() const {return self_twist_flag;}
  INT self_twist_discriminant() const {return CMD;}
  // return base-change code (+1 for base-change, -1 for twisted bc, 0 for neither, 2 for don't know)
  int base_change_code(void) const;
  int cm_code() const {return cm;}
  //  ZZ basis_factor() const {return F->Bdet;}

  // Compute aPmap for first ntp primes if empty, and return it
  map<Quadprime, Eigenvalue> TP_eigs(int ntp, int verbose=0);
  // Compute eQmap if empty and return it, first computing aPmap for
  // first ntp primes if necessary
  map<Quadprime, int> AL_eigs(int ntp=10, int verbose=0);
  // Compute eigmap (principal eigs) if empty and return it
  map<pair<int,string>, Eigenvalue> principal_eigs(int nap, int verbose=0);
  // return the list of traces
  vector<ZZ> traces() const {return trace_list;}

  // NB We only implement file output for newforms with trivial character
  // Filename for this Newform (or conjugate):
  string filename(int conj=0) const;
  // Output newform data (or data for the conjugate newform):
  void output_to_file(int conj=0) const;
  // Input newform data. Returns 0 if data not available, else 1.
  int input_from_file(int verb=0);

  // Construct another newform which is the unramified quadratic twist
  // of this one by D, where D is a discriminant divisor
  Newform unram_quadratic_twist(const INT& D) const;
  // Construct all newforms which are nontrivial unramified quadratic
  // twists of this one, up to Galois conjugacy
  vector<Newform> all_unram_quadratic_twists() const;
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

    // Sort by lists of traces; if these have not been computed yet,
    // or if the character is nontrivial, the trace lists will be
    // empty and we will fall through to the final tests.
    vector<ZZ> traces1 = f1.traces(), traces2 = f2.traces();
    t = traces1<traces2; // true e.g. if f1 has smaller absolute Hecke field degree
    if(t) return 1;
    t = traces1>traces2; // true e.g. if f1 has larger absolute Hecke field degree
    if(t) return 0;

    // We only get here when the traces have not been computed.

    // The characters are the same, we sort by (homological)
    // dimension. NB this is the degree of the principal Hecke field,
    // not the degree of the full Hecke field.
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
  string level_label; // the level's label
  vector<Qideal> Ndivs; // divisors of N
  vector<Quadprime> badprimes; // prime divisrs of N
  vector<Qideal> t2ideals; // list of ideals coprime to level generating 2-torsion in class group
  vector<matop> eps_ops; // list of T(A,A) operators for A in t2ideals

  homspace* H1;  // the ambient modular symbol space at level N
  int cdimH, dimH;
  scalar hmod, dH;

  mat_m T_mat;  // matrix of splitting operator
  string T_name;  // name of splitting operator
  vector<ZZX> factors; // list of multiplicity-1 irreducible factor of charpoly(T)
  // list of possible self-twist discriminants depending only on the level
  vector<INT> possible_self_twists;

  // Return the char poly of T on the new cuspidal subspace using the
  // oldspaces to obtain the old factors with correct multiplicities.

  // If triv_char=0: requires oldspace data for forms with all genus
  // characters, so will only work over fields where this is
  // implemented.

  // If triv_char=1: only requires oldspace data for forms with
  // trivial genus characters, so works over all fields.
  ZZX new_cuspidal_poly(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                        const gmatop &T, int triv_char);


  // Return true iff this combo of ops T has char poly on the new
  // cuspidal subspace which is squarefree and coprime to both the old
  // cuspidal poly and the full Eisenstein poly. f_new is set to the new
  // cuspidal poly. If triv_char=1 then same for the char poly of T on
  // the new cuspidal trivial-character subspace, in which case f_new
  // must also be coprime to the char poly of T on the new cuspidal
  // nontrivial char subspace.
  int valid_splitting_combo(const vector<Quadprime>& Plist, const vector<scalar>& coeffs,
                            const gmatop &T, int triv_char, ZZX& f_new);

  // Find a linear combination T of up to maxnp operators (T_{A,A} or
  // T_P) with coefficients up to maxc, whose char poly on the
  // newspace is squarefree and coprime to its char poly on the
  // oldspace.  Set split_ok=1 if successful else 0.  If
  // triv_char_only, T must have squarefree char poly on the trivial
  // character newspace and be coprime to both its char poly on the
  // full oldspace and on the nontrivial char part of the newspace.
  void find_T(int maxnp, int maxc, int triv_char_only);

public:
  vector<Newform> newforms; // the newforms
  Newspace(void) :verbose(0) {;}
  // constructor from a homspace, looking for a splitting operator
  // using linear combinations of up to maxnp primes, coefficients up
  // to maxc
  Newspace(homspace* h1, int maxnp, int maxc, int triv_char_only, int verb=0);
  int split_ok; // records whether the constructor was able to find a splitting operator

  // constructor from file
  Newspace(const Qideal& level, int verb=0);

  mat_m heckeop(Quadprime& P, int cuspidal=0, int dual=0); // not const as may add info into N
  mat_m heckeop(const matop& T, int cuspidal=0, int dual=0) const;
  mat_m heckeop(const gmatop& T, int cuspidal=0, int dual=0) const;

  int ok() const {return split_ok;}
  int nforms() const {return newforms.size();}
  Qideal level() const {return N;}
  string short_label();
  string long_label();
  string conj_label() const;
  string long_conj_label() const;

  string splitopname() const {return T_name;}

  // Return a list of the degrees of the principal or full Hecke fields
  vector<int> dimensions(int full=0) const;

  // output all newforms: Dimension, Character, Hecke field; optionally aP and AL data
  void display_newforms(int aP=0, int AL=0, int principal_eigs=0, int triv_char_only=0) const;
  // sort the list of newforms using newform_cmp
  void sort_newforms();
  // return the list of newforms
  vector<Newform> the_newforms() const {return newforms;}

  // filename for Newspace (or conjugate)
  string filename(int conj=0);
  // output data for this Newspace (or conjugate) and each Newform
  void output_to_file(int conj=0);
  // Input Newspace data and newform data for each newform. Returns 0 if data missing, else 1.
  int input_from_file(const Qideal& level, int verb=0);

  // For each newform f, create and append all its unramified
  // quadratic twists and resort.  This does nothing if the class
  // number is odd.  We only add twists which are not Galois
  // conjugates: that is, we twist by D in nontrivial cosets of a
  // subgroup of the group of all discriminant divisors.
  void add_unram_quadratic_twists();
};

// dict of Newspaces read from file
extern map<string,Newspace*> Newspace_dict;  // Key: ideal_label(N)
Newspace* get_Newspace(const Qideal& N, int verb=0);

#endif
