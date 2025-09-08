// File NEWFORMS.H

#if     !defined(_NEWFORMS_H)
#define _NEWFORMS_H      1       //flags that this file has been included

#include "ratquads.h"
#include "oldforms.h"
#include "homspace.h"

class newforms;

/* Data stored in a newform data files:
   (Numbers refer to lines of data file)
The first two lines relate to the level:

1. n1ds : number of rational newforms
   n2ds : total dimension of non-rational newforms
   (so n1ds+n2ds is the total new cuspidal dimension)

In even class number only, also:

   new1dims : a list of 2^n2r dimensions summing to n1ds
   new2dims : a list of 2^n2r dimensions summing to n2ds

2. nap  : number of Fourier coefficients

The next 13 lines (3--15) have n1ds (or 2*n1ds) entries each, one per
newform for integers and 2 per newform for quads.  Here, "newform"
really refers to a set of unramified quadratic twists, of size 2^r
(where r=Quad::class_group_2_rank), so just 1 for odd class number, or
for unramified self-twist forms only 2^{r-1}.  The latter condition is
stored in line 16, which is only present when the class number is even.

3.  sfe : sign of functional equation (= - product of aq)
4.  pdot: projection of Manin vector
5.  dp0 : np0=1+Norm(p0)-a_p0, p0 = first good prime
6.  cuspidalfactor : ratio of period / cuspidal period
7.  lambda : twisting prime (or 1), a quad (generator of a principal prime)
8.  lambdadot : scaling factor for lambda
9-12. a, b, c, d : entries of a matrix M=[a,b;N*c,d] in Gamma_0(N) s.t.
13. matdot       : the integral of f over {0,M(0)} is matdot*x
14. bc           : base-change code (see below for definition)
15. cm           : CM discriminant (square-free part), or 0

The next line is only present over field with even class number.
16. CMD: 0, or for forms which are self-twist by an unramified quadratic character, its discriminant.

The next lines contain the Atkin-Lehner eigenvalues for the bad primes
(in standard order):

aq : list of Wq-eigenvalues at bad primes

The last nap lines contain the Fourier Coefficients indexed by primes,
in standard order, including bad primes:

ap : list of Fourier coefficients at all primes

*/

class newform {
friend class newforms;
public:
  newforms *nf;  // pointer to the "parent"
  vec basis;
  vector<long> eigs;   // list of eigenvalues which split off this 1D subspace
  vector<long> aplist; // list of Fourier coefficients, all primes in standard order
  vector<long> aqlist; // list of W-eigenvalues, bad primes in standard order
  scalar dp0;               // 1+N(p0)-a_p0, newforms::p0 = small good principal prime
  scalar pdot;              // Manin vector's projection factor
  rational loverp;       // = pdot/(dp0*nunuits)
  int sfe;               // sign of F.E.
  Quad lambda; scalar lambdadot;  // twisting prime and factor
  Quad a,b,c,d; scalar matdot;    // integration matrix and factor
  int index;             // the index of this newform (from 1)
  int j0; modsym m0; scalar fac, facinv;
  scalar cuspidalfactor;
  INT CMD;            // =D if this is self-twist by unramified disc D dividing Quad::disc, else 0
  vector<long> genus_classes;        // list of genus classes for which we have a nonzero aP
  vector<Qideal> genus_class_ideals; // list of good primes, one in each of these classes
  vector<long> genus_class_aP;       // list of aP for these primes
  int genus_classes_filled;          // set to 1 when genus_classes.size() reaches nchi (when CMD==0)
                                     //or nchi/2 for self-twist (CMD!=0)
  // base-change code:
  //  +1 for b.c. of form over Q with coeffs in Q
  //  +d (d square-free) for b.c. of form over Q with coeffs in Q(sqrt(d))
  //  -1 for twist of b.c. of form over Q with coeffs in Q
  //  -d (d square-free) for twist of b.c. of form over Q with coeffs in Q(sqrt(d))
  //  0  for not b.c. even up to twist
  //  4  for not set
  int bc;
  // CM code:
  //  -d (negative discriminant) if CM with discriminant -d
  //  0 if not CM
  //  1 if not set
  int cm;
  int fake; // flag for "fake rational" forms

  // For each genus class c we count how many primes P in class c have
  // a(P)=0, to aid in detecting self-twist forms.
  vector<int> genus_class_trivial_counter;
  vector<INT> possible_self_twists;

  newform(void) :basis(0), aplist(0) {;}
  // constructor to use just after finding the eigenspace: just sets
  // eigs and basis:
  newform(newforms* nfs, const vec& v, const vector<long>& eigs);

  // constructor to use when data read from file:
  newform(newforms* nfs, int ind,
          const vector<int>& intdata, const vector<Quad>& Quaddata,
          const vector<long>& aq, const vector<long>& ap);

  // After finding all newforms using the basic constructor and
  // setting h1's projcoord and bigtkernbas, fill in the rest of the
  // data for this, extracting the aqlist from eigs:

  //  - aq and sfe from eigs // using compute_AL()
  //  - cuspidalfactor
  //  - L/P                  // using compute_loverp()
  //  - integration matrix   // using find_matrix()

  // Compute AL eigs and SFE
  void compute_AL();
  // Compute cuspidalfactor (needs long basis and bigtkernbas)
  void compute_cuspidalfactor();
  // Compute L/P ratio (needs cuspidalfactor)
  void compute_loverp();
  // Find matrix for integration (needs cuspidalfactor)
  void find_matrix();


  // When a newform has been read from file, we have the aqlist and
  // aplist but not the sequence of eigs in order.  This is needed
  // both for recovering the basis vector from the h1 (in case we want
  // to compute more ap), and for computing oldform multiplcities.
  void eigs_from_data();

  // When a newform has been read from file, when the class number is
  // even,before computing more ap, we need to fill in the genus class
  // data for each newform.
  void fill_in_genus_class_data();

  // For M a *multiple* of this level N, make the list of eigs
  // appropriate for the higher level, deleting the a_P for P dividing
  // M but not N from the sublist of T(P) eigenvalues.
  vector<long> oldform_eigs(Qideal& M);

  // compute the eigenvalue for a single operator on this newform
  // check that the result is factor*a for some a between the bounds
  long eigenvalue(const matop& op, pair<long,long> apbounds, long factor=1);
  // compute aP for this newform, good P
  long eigenvalueHecke(Quadprime& P, int verbose=0);
  // compute A-L eigenvalue for this newform, for the bad prime Q
  long eigenvalueAtkinLehner(Quadprime& Q, int verbose=0);

  void display(void) const;
  void list(string prefix, long nap=-1);

  // Test if form is base-change
  int is_base_change(void);
  // Test if form is base-change up to twist
  int is_base_change_twist(void);
  // if form is base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
  int base_change_discriminant(void);
  // if form is twist of base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
  int base_change_twist_discriminant(void);
  // code to indicate base change or b.c. up to twist
  int base_change_code(void);
  // Test if form is CM, return 0 or the CM disc
  int is_CM(void);
  // Return this twisted by the genus character associated to D
  newform twist(const INT& D);

  // Conjugate data: NB this *only* conjugates data needed for
  // output_to_file(), not everything!
  void conjugate(int debug=0);
};

class newforms :public splitter_base<scalar> {
friend class newform;
private:
  // instantiations of virtual functions required by the splitter_base class:
  mat opmat(int i, int d, int);
  vec opmat_col(int i, int j, int);
  mat opmat_cols(int i, const vec_i& jlist, int);
  mat opmat_restricted(int i, const subspace& s, int d, int);
  smat s_opmat(int i, int d, int);
  svec s_opmat_col(int i, int j, int);
  smat s_opmat_cols(int i, const vec_i& jlist, int);
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int);
  long matdim(void)
  {return h1->dimension;}
  scalar matden(void)
  {return h1->denom3;}

  // For the automatic finding of 1-dimensional egenspaces we need to
  // interface with eclib's splitter class, which wants to know thw
  // i'th operator matrix and its possible eigenvalues, for i=0,1,2,...

  // the list of matrices defining the i'th operator:
  matop h1matop(int);
  // cached list of previous
  vector<matop> h1matops;
  // the list of possible (integer) eigenvalues for the i'th operator:
  vector<long> eigrange(int i);
  // cached list of previous
  vector<vector<long>> eigranges;

  long dimoldpart(const vector<long> l) {return of->dimoldpart(l);}

  // Compute a long basis vector from a short one
  vec lengthen_basis(const vec& sbasis);

  long j0; // single j0 in 1..ngens, a pivot for all newforms, or 0
  std::set<long> jlist; // set of j in 1..ngens, including pivots for all newforms
  map<long,modsym> mjlist; // corresponding modular symbols
// Look for a pivotal index j0 (from 1) such that
  // nflist[i].basis[j0]!=0 for all i, or a set of such j (Each
  // newform stores a j0-value and this nonzero coordinate as "fac".)
  void find_jlist();

protected:
  oldforms *of;      // pointer to one, not an array
  Quadprime P0;      // the first good principal prime
  int iP0;           // P0's index in the list goodprimes (indexed from 0)
  scalar nP0;        // N(P0)
  vector<long> aP0;  // vector of ap for P0 (indexed from 0)
  vec mvp;           // Manin vector of P0
  vec zero_infinity; // {0,oo}
public:
  Qideal N;  // the level
  vector<Quadprime> badprimes; // list of all bad primes Q (dividing the level N)
  vector<Quadprime> goodprimes;  // good primes in order
  vector<int> bad_prime_conjugation_permutation;
  vector<Qideal> nulist; // list of ideals coprime to level generating 2-torsion in class group
  vector<INT> possible_self_twists;
  int level_is_square;
  int verbose, nwq, nap, n2r, nchi;

  // subdimensions of rational/non-rational new trivial character cuspidal subspaces
  int n1ds, n2ds;
  // dimensions of the trivial character cuspidal subspace, with its old and new subspaces:
  int dimtrivcusp, dimtrivcuspold, dimtrivcuspnew;
  // the same divided into 2^n2r parts, the i'th part being the
  // dimension of the subspace which is self-twist by the i'th
  // unramified quadratic character, except for i=0 for the part which
  // is not self-twist.  (Ususally all the dimension is in this -'th
  // part, but not always).
  vector<int> alldims, olddims, newdims;
  // Each of the previous dimensions is divided into a 1-dimensional
  // part (spanned by rational forms) and a >=2-dimensional part (the
  // rest). In particular, n1ds=sum(new1dims) and n2ds=sum(new2dims).
  // Ususally these dims are all concentrated in the 0'th component
  // (and always in odd class number when this is the only component).
  vector<int> old1dims, new1dims;
  vector<int> old2dims, new2dims;
  homspace* h1; // pointer to one, not an array
  scalar hmod, modulus, characteristic; // 0 or prime
  int have_bases;
  vector<newform> nflist;
  explicit newforms(const Qideal& N, scalar mod, int disp=0, scalar ch=scalar(0));
  ~newforms(void)
  {
    if(h1!=NULL)delete h1;
    if(of!=NULL)delete of;
  }
  void display(int detail=1);
  // List newforms in a fixed format. NB In even class number, each
  // element from nflist gives rise to 2**n2r (or 2**(n2r-1) in case
  // of forms with self-twist) actual rational newforms.
  void list(long nap=-1);

private:
  // Compute the associated homspace
  void makeh1(void);
  // Set projcoord member of homspace
  void make_projcoord(void);
  // Set bigtkernbas member of homspace
  void make_bigtkernbas(void);
  // fill in extra data for all:
  // projcoord, jlist, zero_infinity, mvp, lambdas
  // and, for each newform:
  // AL-eigenvalues and SFE (if AL);
  // loverp if LP;
  // integration matrix and multiple (if M);
  void fill_in_newform_data(int AL=1, int CF=1, int LP=1, int M=1);
  void find_lambdas();

  // add newform with basis b1, eiglist eigs to current list (b2 not used)
  // if use_nf_number=-1;  if >=0 just store eigs and basis in that preexisting nf
  int use_nf_number;
  void use(const vec& b1, const vec& b2, const vector<long> eigs);

 public:

  // methods for computing Hecke eigenvalues

  void getap(int first, int last, int verbose=0);

  // compute list of eigenvalues at P for each newform
  vector<long> apvec(Quadprime& P);
  // compute list of eigenvalues at P for each newform (good P, Euclidean) and check that it is within bounds
  vector<long> apvec_euclidean(Quadprime& P, pair<long,long> apbounds);
  // compute list of eigenvalues of op for each newform and check that it is within bounds
  vector<long> apvec(const matop& op, pair<long,long> apbounds);
  // compute list of eigenvalues given the image images[j] for each j in jlist
  vector<long> apvec_from_images(map<int,vec> images, pair<long,long> apbounds, const string& name);

  void output_to_file(string eigfile) const;

  // sorting functions
  void sort_eigs(void);
  void sort_lmfdb(void);

  // find newforms by splitting homspace
  void find();
  // try to read from file, return 1 for success (data file exists) or 0 if no data file exists
  int read_from_file();
  // try to read from file, and if no data file exists, finds from scratch and stores
  void read_from_file_or_find();
  // if created from stored data but need bases and homspace
  // (if extra_data, also recompute cuspidalfactor, loverp, integration data)
  void makebases(int extra_data=1);

  // Conjugate data: NB this *only* conjugates data needed for
  // output_to_file(), not everything!
  void conjugate(int debug=0);
};

#endif
