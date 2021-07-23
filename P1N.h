// FILE P1N.H: declaration of class for P^1(O_K/N) for an arbitrary ideal N

#if     !defined(_P1N_H)
#define _P1N_H      1       //flags that this file has been included

#include <iostream>

#include <eclib/arith.h>
#include "qideal.h"
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
  P1N(const Qideal& I);                                //constructor
  void make_symb(long i, Quad& c, Quad& d); // assign c, d to the i'th (c:d) symbol
  long index(const Quad& c, const Quad& d); // index i of (c:d)
  long size() {return psi;}

  long merge_indices(const vector<long>& klist) {return ::merge_indices(psilist, klist);}
  vector<long> split_indices(long k) {return ::split_indices(psilist, k);}

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


#endif
