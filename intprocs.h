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

// define INT_IS_ZZ or INT_IS_long to use ZZ or long int as base integer type for Quads

#ifdef INT_IS_long
typedef long INT; // integer type for components of a Quad
#endif
#ifdef INT_IS_ZZ
#include <eclib/marith.h>
typedef bigint INT; // integer type for components of a Quad
inline bigfloat to_bigfloat(const bigint& x) { return to_RR(x);}
INT rounded_division(INT aa, INT bb);
int divides(const INT& aa, const INT& bb);
#endif

const INT ZERO(0);
const INT ONE(1);
const INT TWO(2);
const INT THREE(3);

inline long I2long(long n) {return n;}
void sqrt_mod_p(long & x, long a, long p);
int divides(long a, long b, long& q, long& r);
long squarefree_part(long d);

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
//
long rounded_division(long aa, long bb);

//functions needed for non-euclidean fields to compute bezout/quadgcd

// content
INT vecgcd(const vector<INT>& a);

//returns g = content(a) = a.c
INT vecbezout(const vector<INT>& a, vector<INT>& c);

// dot product
INT dot(const vector<INT>& a, const vector<INT>& c);

// Finds basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

void findzbasis(const vector<INT>& first, const vector<INT>& second, vector<INT>& basis);

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
inline long from_bits(vector<int> aa) {return from_bits(aa, aa.size());}

// For D1, D fundamental discriminants, test if D=D1*D2 with D2 another discriminant
int div_disc(INT D1, INT D);

// return list of integers from first to last inclusive
vector<long> range(long first, long last);

// return 1 with r=sqrt(a) if a is square, else return 0:
long is_square(long a, long& r);

vec reduce_modp(const vec& v, const scalar& p=DEFAULT_MODULUS);
mat reduce_modp(const mat& m, const scalar& p=DEFAULT_MODULUS);

#endif
