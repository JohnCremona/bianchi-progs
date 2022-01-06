// FILE P1N.H: declaration of class for P^1(O_K/N) for an arbitrary ideal N

#if     !defined(_P1N_H)
#define _P1N_H      1       //flags that this file has been included

#include <iostream>

#include "mat22.h"
#include "primes.h"

// utilities for a standard bijection between [0,1,...,n-1] and the
// product of [0,1,...,n_i-1] where n is the product of the n_i

long merge_indices(const vector<long>& nlist, const vector<long>& klist);
vector<long> split_indices(const vector<long>& nlist, long k);

class P1N {
  vector<long> residue_codes, noninvertible_residue_indices;

  // for prime powers only (NB resnum[0]=0 always):
  // : residue_codes[i] = j >0 if resnum[i] is invertible with inverse resnum[j]
  // : residue_codes[i] =-j <=0 if resnum[i] is the j'th noninvertible residue
  // : noninvertible_residue_indices[i] = j if the i'th noninvertible residue is resnum[j]

  vector<P1N> P1PP; // one for each prime power dividing N when N is *not* a prime power
  vector<long> psilist; // phi of each prime power
public:
  P1N() {;}                                            //constructor
  explicit P1N(const Qideal& I);                                //constructor

  P1N& operator=(const P1N& other)
  {
    residue_codes = other.residue_codes;
    noninvertible_residue_indices = other.noninvertible_residue_indices;
    P1PP = other.P1PP;
    psilist = other.psilist;
    nrm = other.nrm;
    phi = other.phi;
    psi = other.psi;
    np = other.np;
    N = other.N;
    return *this;
  }
  void make_symb(long i, Quad& c, Quad& d); // assign c, d to the i'th (c:d) symbol
  void reduce(Quad& c, Quad& d);            // simplify c,d without changing (c:d)
  long index(const Quad& c, const Quad& d); // index i of (c:d)
  long size() const {return psi;}
  Qideal level() const {return N;}

  long merge_indices(const vector<long>& klist) const {return ::merge_indices(psilist, klist);}
  vector<long> split_indices(long k) const {return ::split_indices(psilist, k);}

  // return a matrix M = [a, b; c, d] with det=1 lifting the i'th (c:d) symbol
  mat22 lift_to_SL2(long i);

  // each M in GL2 permutes the (c:d) symbols by right multiplcation:
  long apply(const mat22& M, long i);

  // test function
  void check(int verbose=0);
protected:
  long nrm, phi, psi;
  Qideal N;  // the level
  int np;    // number of bad primes
};

// compute a matrix M = [a, b; c, d] with det=1 lifting (c:d) in P^1(N)
mat22 lift_to_SL2(Qideal& N, const Quad& cc, const Quad& dd);

// class for action of 2x2 matrices on a P1N
class action  :public mat22 {
private:
  P1N* level;
public:
  action() {;}
  action(P1N* N, const mat22& mm) : mat22(mm), level(N) {}
  action(P1N* N, const Quad& a, const Quad& b, const Quad& c, const Quad& d) : mat22(a,b,c,d), level(N)  {}

  // Right action of a 2x2 matrix on P^1(N) mapping input symbol's
  // index to output symbol's index:

  int operator()(int i) const
  {
    int j = level->apply(*this, i);
    //    cout<<" - applying "<<(*this)<<" to symbol #"<<i<<" gives symbol #"<<j<<endl;
    return j;
  }
};

#endif
