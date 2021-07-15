// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include <eclib/arith.h>
// For b>0, roundover(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
long roundover(long aa, long bb);

// eclib only has the following for bigints
long sqrt_mod_p(long a, long p); // p odd prime, a quadratic residue

//functions needed for non-euclidean fields to compute bezout/quadgcd

// content
long vecgcd(const vector<long>& a);

//returns g = content(a) = a.c
long vecbezout(const vector<long>& a, vector<long>& c);

// dot product
long dot(const vector<long>& a, const vector<long>& c);

// sets basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis
// for the Z-module spanned by [first[i], second[i]], and also sets x, y to be vectors such that
//
// first.x = e1
// first.y = e2
// second.x = f1
// second.y = 0

void findzbasiscoeffs(const vector<long>& first, const vector<long>& second,
                      vector<long>& basis, vector<long>& x, vector<long>& y);

//Same as findzbasiscoeffs except don't need x,y: sets
//basis={e1,e2,f1} such that [[e1,f1], [e2,0]] is a Z-basis for the
//Z-module spanned by [first[i], second[i]]

void findzbasis(const vector<long>& first, const vector<long>& second, vector<long>& basis);

#endif
