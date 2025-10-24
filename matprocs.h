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
template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A);

// convert an eclib mat to an NTL mat_ZZ_p:
mat_ZZ_p mat_to_mat_ZZ_p(mat A);

// return f(X/c)*c^d
ZZX scale_poly_up(const ZZX& f, const ZZ& c);
// return f(c*X)/c^d
ZZX scale_poly_down(const ZZX& f, const ZZ& c);

// return f(X) mod m (or just f if m==0)
ZZX reduce_poly(const ZZX& f, const ZZ& m);

// Coprime test
int AreCoprime(const ZZX& f, const ZZX& g);

// Monic test
inline int IsMonic(const ZZX& f)
{ return IsOne(LeadCoeff(f));}

// Irreducibility test (ignoring content)
int IsIrreducible(const ZZX& f);

// Squarefree test
int IsSquareFree(const ZZX& f);

// compute char poly of A:
ZZX charpoly(const mat_ZZ& A);

// compute char poly of A:
inline ZZX charpoly(const mat_m& A) {return charpoly(mat_to_mat_ZZ(A));}

// compute char poly of A/den mod m:
ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den = to_ZZ(1), const scalar& m = scalar(0));

// return f(X^2)
ZZX XtoX2(const ZZX& f);

// write f(x) = f0(X^2)+X*f1(X^2)
void parity_split(const ZZX& f, ZZX& f0, ZZX& f1);

// assuming f irreducible:
// return 0 if f(x^2) is irreducible; else
// return 1 and set g where f(x^2)=g(x)g(-x) (*-1 if degree odd)
int is_square(const ZZX& f, ZZX& g);

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m);

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A);
// evaluate f(A) (assumes f monic)
template<class T>
mat_m evaluate(const ZZX& f, const Zmat<T>& A);

// p should be monic:
mat_ZZ CompanionMatrix(const ZZX& p);

// check that a matrix is a scaled involution (modulo m, if m!=0):
int check_involution(const mat_ZZ& A, scalar den, const scalar& m, int verbose=0);

// check that two matrices commute (modulo m, if m!=0):
int commute(const mat_ZZ& A, const mat_ZZ& B, const scalar& m);

// check that a matrix commutes (modulo m, if m!=0) with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const scalar& m);

// factor a primitive (e.g. monic) polynomial
vec_pair_ZZX_long factor(const ZZX& f);

// display factors of a polynomial:
void display_factors(const ZZX& f);

// display factors of a polynomial mod p:
void display_factors(const ZZ_pX& f);

// rank of an NTL matrix:
long rank(mat_ZZ A);

// nullity of an NTL matrix:
long nullity(mat_ZZ A);

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats);
// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats);

// same as m.output(cout) except no newlines between rows
template<class T>
void output_flat_matrix(const Zmat<T>& m, ostream&s = cout);

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
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(),
                                        fac2.a.rep.begin(), fac2.a.rep.end());
  }
};

//bool operator<(ZZ_p a, ZZ_p b); // {return rep(a)<rep(b);}

struct factor_modp_comparison {
  bool operator()(pair_ZZ_pX_long& fac1, pair_ZZ_pX_long& fac2)
  {
    auto cmp = [](const ZZ_p& a, const ZZ_p& b) {return rep(a)<rep(b);};
    // first sort by degree of the factor
    int s = deg(fac1.a) - deg(fac2.a);
    if(s) return (s<0); // true if fac1 has smaller degree

    // then sort by exponent of the factor
    s = fac1.b - fac2.b;
    if(s) return (s<0); // true if fac1 is to a lower exponent

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(fac1.a.rep.begin(), fac1.a.rep.end(),
                                        fac2.a.rep.begin(), fac2.a.rep.end(),
                                        cmp);
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
    return std::lexicographical_compare(pol1.rep.begin(), pol1.rep.end(),
                                        pol2.rep.begin(), pol2.rep.end());
  }
};

struct poly_modp_comparison {
  bool operator()(ZZ_pX& pol1, ZZ_pX& pol2)
  {
    auto cmp = [](const ZZ_p& a, const ZZ_p& b) {return rep(a)<rep(b);};
    // first sort by degree of the factor
    int s = deg(pol1) - deg(pol2);
    if(s) return (s<0); // true if pol1 has smaller degree

    // finally lexicographically compare the coefficient lists
    return std::lexicographical_compare(pol1.rep.begin(), pol1.rep.end(),
                                        pol2.rep.begin(), pol2.rep.end(),
                                        cmp);
  }
};

extern factor_comparison fact_cmp;
extern poly_comparison poly_cmp;
extern factor_modp_comparison fact_modp_cmp;
extern poly_modp_comparison poly_modp_cmp;

vector<ZZ> coeffs(const ZZX& p);
vector<ZZ> coeffs(const ZZ_pX& p);
string polynomial_string(const vector<ZZ>& coeffs, const string& var="X");
string polynomial_string(const vec_m& coeffs, const string& var="X");
string polynomial_string(const ZZX& p, const string& var="X");
string polynomial_string(const ZZ_pX& p, const string& var="X");

#endif
