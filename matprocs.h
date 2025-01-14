// FILE MATPROCS.H:  interface with NTL mat_ZZ and related classes

#if     !defined(_MATPROCS_H)
#define _MATPROCS_H      1       //flags that this file has been included

// include stuff from eclib
#include "intprocs.h"
#include <eclib/matrix.h>

// include stuff from NTL
#include <NTL/LLL.h>

#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>

#include <NTL/mat_poly_ZZ_p.h>
#include <NTL/ZZ_pXFactoring.h>

#include "qidloop.h"

// convert an eclib mat to an NTL mat_ZZ:
mat_ZZ mat_to_mat_ZZ(mat A);

// convert an eclib mat to an NTL mat_ZZ_p:
mat_ZZ_p mat_to_mat_ZZ_p(mat A);

// compute char poly of A/den:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den);

// check that a matrix is a scaled involution:
int check_involution(const mat_ZZ& A, long den=1, int verbose=0);

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist);

// display factors of a polynomial:
void display_factors(const ZZX& f);

// rank of an NTL matrix:
long rank(mat_ZZ A);

// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// function to sort a factorization vector, first by degree of factor
// then exponent of factor then lexicographically

struct factor_comparison {
  bool operator()(pair_ZZX_long& fac1, pair_ZZX_long& fac2)
  {
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(), fac2.a.rep.begin(), fac2.a.rep.end());
  }
};

extern factor_comparison fact_cmp;

#endif
