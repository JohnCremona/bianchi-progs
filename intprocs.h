// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include "arith_extras.h"

// EITHER define INT_IS_ZZ or INT_IS_long to use ZZ or long int as base integer type for Quads
// OR define FLINT to use FLINT's fmpz

#ifdef FLINT
#include "int.h"
#endif

#ifdef INT_IS_long
#include <eclib/marith.h>
typedef long INT; // integer type for components of a Quad
#endif

#ifdef INT_IS_ZZ
#include <eclib/marith.h>
typedef bigint INT; // integer type for components of a Quad
inline bigfloat to_bigfloat(const bigint& x) { return to_RR(x);}
// Return q=round(a/b) with halves going down if round_down=1 (default) else up
INT rounded_division(INT aa, INT bb, int round_down=1);
int divides(const INT& aa, const INT& bb);
int divrem(const bigint& a, const bigint& b, bigint& q, bigint& r);
INT isqrt(const INT& a);
#endif

const INT ZERO(0);
const INT ONE(1);
const INT MONE(-1);
const INT TWO(2);
const INT THREE(3);

// content
INT content(const vector<INT>& a);

//returns g = content(a) = a.c
INT vecbezout(const vector<INT>& a, vector<INT>& c);

// dot product
INT dot(vector<INT>& a, vector<INT>& c);

// Test if the Z-module spanned by [coords[0][i], coords[1][i]] is all of Z^2 (for coprimality testing)
int span_Z2( const pair< vector<INT>, vector<INT>> & coords);

// Finds basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

vector<INT> Zbasis(const pair<vector<INT>, vector<INT>>& coords);

// For D1, D fundamental discriminants, test if D=D1*D2 with D2 another discriminant
int div_disc(INT D1, INT D);

#endif
