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
  mat opmat(int i, int d, int v=0)
  {return h1->opmat(i,d,v);}
  vec opmat_col(int i, int j, int v=0)
  {return h1->opmat_col(i,j,v);}
  mat opmat_cols(int i, const vec& jlist, int v=0)
  {return h1->opmat_cols(i,jlist,v);}
  mat opmat_restricted(int i, const subspace& s, int d, int v=0)
  {return h1->opmat_restricted(i,s,d,v);}
  smat s_opmat(int i, int d, int v=0)
  {return h1->s_opmat(i,d,v);}
  svec s_opmat_col(int i, int j, int v=0)
  {return h1->s_opmat_col(i,j,v);}
  smat s_opmat_cols(int i, const vec& jlist, int v=0)
  {return h1->s_opmat_cols(i,jlist,v);}
  smat s_opmat_restricted(int i, const ssubspace& s, int d, int v=0)
  {return h1->s_opmat_restricted(i,s,d,v);}
  long matdim(void);
  long matden(void);
  vector<long> eigrange(int i) {return h1->eigrange(i);}
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
  vector<Quadprime> plist; // bad primes
  int is_square, npdivs;
  int verbose, n1ds,n2ds, nnflist, nap, ntp, nwq;
  homspace* h1; // pointer to one, not an array
  long hmod, nfhmod;
  vector<newform> nflist;
  explicit newforms(const Qideal& N, int disp=0);
  ~newforms(void) {
                   if(h1)delete h1;
                  }
  void init();
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
  vector<long> apvec(Quadprime& P);  // computes a[P] for each newform

  void output_to_file(string eigfile) const;

  // sorting functions
  void sort_eigs(void);
  void sort_lmfdb(void);

  void createfromscratch();
  void createfromdata();
  void makebases(); // if created from stored data but need bases and homspace
};

#endif
