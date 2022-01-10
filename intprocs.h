// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include <eclib/interface.h>
#include <eclib/marith.h>

//#define QUINT_IS_ZZ

#ifndef QUINT_IS_ZZ
#define QUINT_IS_long
typedef long QUINT; // integer type for components of a Quad
#else
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

#endif
