// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

// for convenience we put all the includes from eclib here, since all
// other files include this one directly or indirectly:

#include <eclib/xsplit.h>   // which includes method.h
#include <eclib/rat.h>
#include <eclib/unimod.h>
#include <eclib/timer.h>
#include <eclib/curvesort.h> // for letter codes

// define this in the Makefile (and make clean) to use ZZ as base integer type instead of long int
//#define QUINT_IS_ZZ

#ifndef QUINT_IS_ZZ
#define QUINT_IS_long
typedef long QUINT; // integer type for components of a Quad
#else
#include <eclib/marith.h>
typedef bigint QUINT; // integer type for components of a Quad
inline bigfloat to_bigfloat(const bigint& x) { return to_RR(x);}
#endif

const QUINT ZERO(0);
const QUINT ONE(1);
const QUINT TWO(2);
const QUINT THREE(3);

inline long I2long(long n) {return n;}
inline int div(long a, long b) {return (a==0? b==0: (b%a)==0);}
inline int is_nonnegative(long a) {return a>=0;}
inline int is_positive(long a) {return a>0;}
void sqrt_mod_p(long & x, long a, long p);
int divides(long a, long b, long& q, long& r);
long squarefree_part(long d);

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
//
#ifdef QUINT_IS_ZZ
int is_nonnegative(QUINT a);
QUINT rounded_division(QUINT aa, QUINT bb);
#endif
long rounded_division(long aa, long bb);

//functions needed for non-euclidean fields to compute bezout/quadgcd

// content
QUINT vecgcd(const vector<QUINT>& a);

//returns g = content(a) = a.c
QUINT vecbezout(const vector<QUINT>& a, vector<QUINT>& c);

// dot product
QUINT dot(const vector<QUINT>& a, const vector<QUINT>& c);

// Finds basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

void findzbasis(const vector<QUINT>& first, const vector<QUINT>& second, vector<QUINT>& basis);

// return the i'th bit of a (a>=0) (not needed, NTL has this)
// inline int bit(long a, long i) {return (a& (1<<i));}
// return the dot product (0/1) of a and b bitwise, with 0<=a,b<2^r
int dotbits(long a, long b, int r);
// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r);
// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(vector<long> alist, int r);
// return list of bits of a
vector<int> bits(long a, int r);
// recover a from its bit vector of length r
long from_bits(vector<int> aa, int r);

// return list of integers from first to last inclusive
vector<long> range(long first, long last);

// return 1 with r=sqrt(a) if a is square, else return 0:
long is_square(long a, long& r);

vec reduce_modp(const vec& v, const scalar& p=DEFAULT_MODULUS);
mat reduce_modp(const mat& m, const scalar& p=DEFAULT_MODULUS);

#endif
