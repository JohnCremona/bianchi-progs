// FILE ARITH_EXTRAS.H: miscellaneous integer (int/long) functions which could be in eclib/arith.h

#if     !defined(_ARITH_EXTRAS_H)
#define _ARITH_EXTRAS_H      1       //flags that this file has been included

// for convenience we put all the includes from eclib here, since all
// other files include this one directly or indirectly:

#include <eclib/xsplit.h>   // which includes method.h
#include <eclib/rat.h>
#include <eclib/unimod.h>
#include <eclib/timer.h>
#include <eclib/curvesort.h> // for letter codes

inline long I2long(long n) {return n;}
int isqrt(long a, long& root);
long isqrt(const long a);
int divrem(long a, long b, long& q, long& r);
long squarefree_part(long d);
void sqrt_mod_p(long & x, long a, long p);

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
//
long rounded_division(long aa, long bb);

// return list of integers from first to last inclusive
vector<long> range(long first, long last);

// return 1 with r=sqrt(a) if a is square, else return 0:
long is_square(long a, long& r);

// return the dot product (0/1) of a and b bitwise, with 0<=a,b<2^r
int dotbits(long a, long b, int r);
// return list of bits of a
vector<int> bits(long a, int r);
// recover a from its bit vector of length r
long from_bits(vector<int> aa, int r);
inline long from_bits(vector<int> aa) {return from_bits(aa, aa.size());}

// return a basis for the orthogonal complement of a<2^r (viewed as a bit vector of length r)
vector<long> dotperp(long a, int r);
// return a basis for the orthogonal complement of the span of a in alist (viewed as bit vectors of length r)
vector<long> dotperp(vector<long> alist, int r);

vec reduce_modp(const vec& v, const scalar& p=DEFAULT_MODULUS);
mat reduce_modp(const mat& m, const scalar& p=DEFAULT_MODULUS);

#endif
