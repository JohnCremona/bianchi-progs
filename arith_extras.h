// FILE ARITH_EXTRAS.H: miscellaneous integer (int/long) functions

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

void sqrt_mod_p(long & x, long a, long p);

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

#endif
