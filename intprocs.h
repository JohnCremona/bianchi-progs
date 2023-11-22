// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include "arith_extras.h"

// define INT_IS_ZZ or INT_IS_long to use ZZ or long int as base integer type for Quads

#ifdef FLINT
#include "flint.h"
#else

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
#endif

const INT ZERO(0);
const INT ONE(1);
const INT MONE(-1);
const INT TWO(2);
const INT THREE(3);

//functions needed for non-euclidean fields to compute bezout/quadgcd

// content
INT vecgcd(const vector<INT>& a);

//returns g = content(a) = a.c
INT vecbezout(vector<INT>& a, vector<INT>& c);

// dot product
INT dot(vector<INT>& a, vector<INT>& c);

// Finds basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

void findzbasis(const vector<INT>& first, const vector<INT>& second, vector<INT>& basis);

// For D1, D fundamental discriminants, test if D=D1*D2 with D2 another discriminant
int div_disc(INT D1, INT D);

#endif
