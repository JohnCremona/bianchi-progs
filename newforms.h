// File NEWFORMS.H

#if     !defined(_NEWFORMS_H)
#define _NEWFORMS_H      1       //flags that this file has been included

#include <eclib/xsplit.h>   // which includes method.h
#include <eclib/rat.h>

#include "ratquads.h"
#include "oldforms.h"
#include "homspace.h"

class newforms;

/* Data stored in a newform data files:
   (Numbers refer to lines of data file)
The first two lines relate to the level:

1. n1ds : number of rational newforms
   n2ds : total dimension of non-rational newforms
2. nap  : number of Fourier coefficients

The next 11 lines (3--13) have n1ds (or 2*n1ds) entries each, one per
newform for integers and 2 per newform for quads:

3.  sfe : sign of functional equation (= - product of aq)
4.  pdot: projection of Manin vector
5.  dp0 : np0=1+Norm(p0)-a_p0, p0 = first good prime
6.  cuspidalfactor : ratio of period / cuspidal period
7.  lambda : twisting prime (or 1), a quad (generator of a principal prime)
8.  lambdadot : scaling factor for lambda
9-12. a, b, c, d : entries of a matrix M=[a,b;N*c,d] in Gamma_0(N) s.t.
13. matdot       : the integral of f over {0,M(0)} is matdot*x

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
  int dp0;               // 1+N(p0)-a_p0, newforms::p0 = small good principal prime
  int pdot;              // Manin vector's projection factor
  rational loverp;       // = pdot/(dp0*nunuits)
  int sfe;               // sign of F.E.
  Quad lambda; int lambdadot;  // twisting prime and factor
  Quad a,b,c,d; int matdot;    // integration matrix and factor
  int j0; int fac, facinv;
  long cuspidalfactor;

  newform(void) :basis(0), aplist(0) {;}
  // constructor to use just after finding the eigenspace: just sets
  // eigs and basis:
  newform(newforms* nfs, const vec& v, const vector<long>& eigs);

  // constructor to use when data read from file:
  newform(newforms* nfs,
          const vector<int>& intdata, const vector<Quad>& Quaddata,
          const vector<long>& aq, const vector<long>& ap);

  // After finding all newforms using the basic constructor and
  // setting h1's projcoord, fill in the rest of the data for this,
  // given that it is the j'th newform:
  // - ap and sfe
  // - L/P
  // - manin vector data
  // - integration matrix using find_matrix()

  void fill_in_data(int j);

  void display(void) const;
  void list(long nap=-1) const;
  // To find cuspidal factor:
  void find_cuspidal_factor(const vec& v);
  // To find matrix for integration:
  void find_matrix(int j);
  // Test if form is base-change
  int is_base_change(void) const;
  // Test if form is base-change up to twist
  int is_base_change_twist(void) const;
  // if form is base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
  int base_change_discriminant(void) const;
  // if form is twist of base-change, find the d s.t. the bc has eigenvalues in Q(sqrt(d))
  int base_change_twist_discriminant(void) const;
  // Test if form is CM, return 0 or the CM disc
  int is_CM(void) const;
};

class newforms :public splitter_base {
friend class newform;
private:
  int dimsplit, maxdepth, upperbound;

  // instantiations of virtual functions required by the splitter_base class:
  mat opmat(int i, int d, int v=0);
  vec opmat_col(int i, int j, int v=0);
  mat opmat_cols(int i, const vec& jlist, int v=0);
  mat opmat_restricted(int i, const subspace& s, int d, int v=0);
  smat s_opmat(int i, int d, int v=0);
  svec s_opmat_col(int i, int j, int v=0);
  smat s_opmat_cols(int i, const vec& jlist, int v=0);
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0);
  long matdim(void)
  {return h1->dimension;}
  long matden(void)
  {return h1->denom3;}

  // For the automatic finding of 1-dimensional egenspaces we need to
  // interface with eclib's splitter class, which wants to know thw
  // i'th operator matix and its possible eigenvalues, for i=0,1,2,...

  // the list of matrices defining the i'th operator:
  matop h1matop(int);
  // the list of possible (integer) eigenvalues for the i'th operator:
  vector<long> eigrange(int i);

  long dimoldpart(const vector<long> l) {return of->dimoldpart(l);}

  // data used for ap computation
  int easy;
  vector<long> pdotlist, pdotlistinv;

  long j0;
  std::set<long> jlist;
  // Look for a pivotal index j0 (from 1) such that
  // nflist[i].basis[j0]!=0 for all i, or a set of such j (Each
  // newform stores a j0-value and this nonzero coordinate as "fac".)
  void find_jlist();

protected:
  oldforms *of; // pointer to one, not an array
  Quadprime P0; int iP0; long nP0; vec mvp;
  vector<long> aP0;
  vec zero_infinity;
public:
  Qideal N;  // the level
  vector<Quadprime> badprimes; // list of bad primes Q with square ideal class or even exponent
  vector<Quadprime> goodprimes;  // good primes in order
  vector<Qideal> nulist; // list of ideals coprime to level generating 2-torsion in class group
  int is_square;
  int verbose, n1ds,n2ds, nnflist, nwq, nap, n2r;
  homspace* h1; // pointer to one, not an array
  long hmod, nfhmod;
  long characteristic; // 0 or prime
  vector<newform> nflist;
  explicit newforms(const Qideal& N, int disp=0, long ch=0);
  ~newforms(void) {
                   if(h1)delete h1;
                  }
  void display(int detail=1);
  void list(long nap=-1);

private:
  // Compute the associated homspace
  void makeh1plus(void);
  // Set projcoord member of homspace
  void make_projcoord(void);
  // fill in extra data in each newforms:
  void fill_in_newform_data(int everything=1);
  void find_lambdas();

  // add newform with basis b1, eiglist eigs to current list (b2 not used)
  // if use_nf_number=-1;  if >=0 just store eigs and basis in that preexisting nf
  int use_nf_number;
  void use(const vec& b1, const vec& b2, const vector<long> eigs);

 public:

  // methods for computing Hecke eigenvalues

  void getap(int first, int last, int verbose=0);
  void getoneap(Quadprime& P, int verbose=0, int store=1);
  // compute eigenvalue at P for each newform
  vector<long> apvec(Quadprime& P);
  // compute eigenvalue at P for each newform (good P, Euclidean) and check that it is in elist
  vector<long> apvec_euclidean(Quadprime& P, const vector<long>& elist);
  // compute eigenvalue of op for each newform and check that it is in elist
  vector<long> apvec(const matop& op, const vector<long>& elist);

  void output_to_file(string eigfile) const;

  // sorting functions
  void sort_eigs(void);
  void sort_lmfdb(void);

  void createfromscratch();
  void createfromdata();
  void makebases(); // if created from stored data but need bases and homspace
};


// List of bad primes (dividing N) followed by good primes to length
// at least np, making sure that the list includes at least one good
// principal prime.  iP0 is set to the index in the list of the first
// good principal prime.
// If p (default 0) is nonzero, omit bad primes and primes dividing P
vector<Quadprime> make_primelist(Qideal& N, int np, int& iP0, int p=0);

// compute a list of ideals coprime to N whose classes generate the 2-torsion
vector<Qideal> make_nulist(Qideal& N);
// compute a list of primes Q dividing N with Q^e||N such that [Q^e] is square
vector<Quadprime> make_badprimes(Qideal& N);
// compute a list of at least nap good primes (excluding those
// dividing characteristic if >0), to include at least on principal
// one which has index iP0;
vector<Quadprime> make_goodprimes(Qideal& N,  int np, int& iP0, int p=0);


#endif
