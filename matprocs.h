// FILE MATPROCS.H:  interface with NTL mat_ZZ and related classes

#if     !defined(_MATPROCS_H)
#define _MATPROCS_H      1       //flags that this file has been included

// include stuff from eclib
#include "intprocs.h"

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

// return f(X/c)*c^d
ZZX scale_poly_up(const ZZX& f, const ZZ& c);
// return f(c*X)/c^d
ZZX scale_poly_down(const ZZX& f, const ZZ& c);

// return f(X) mod m (or just f if m==0)
ZZX reduce_poly(const ZZX& f, const ZZ& m);

// Squarefree test
int IsSquareFree(const ZZX& f);

// compute char poly of A/den mod m:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den, const scalar& m);

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m);

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A);
// evaluate f(A) (assumes f monic)
mat evaluate(const ZZX& f, const mat& A);

// check that a matrix is a scaled involution (modulo m, if m!=0):
int check_involution(const mat_ZZ& A, scalar den, const scalar& m, int verbose=0);

// check that two matrices commute (modulo m, if m!=0):
int commute(const mat_ZZ& A, const mat_ZZ& B, const scalar& m);

// check that a matrix commutes (modulo m, if m!=0) with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const scalar& m);

// display factors of a polynomial:
void display_factors(const ZZX& f);

// rank of an NTL matrix:
long rank(mat_ZZ A);

// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// function to sort factorizations (lists of (factor,exponent) pairs),
// first by degree of factor then exponent of factor then
// lexicographically

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

// function to sort of polynomials, first by degree of factor
// then lexicographically

struct poly_comparison {
  bool operator()(ZZX& pol1, ZZX& pol2)
  {
    // first sort by degree of the factor
    int s = deg(pol1) - deg(pol2);
    if(s) return (s<0); // true if pol1 has smaller degree

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(pol1.rep.begin(), pol1.rep.end(), pol2.rep.begin(), pol2.rep.end());
  }
};

extern factor_comparison fact_cmp;
extern poly_comparison poly_cmp;

#endif
