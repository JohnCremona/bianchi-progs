// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include <eclib/arith.h>

typedef long QUINT; // integer type for components of a Quad

long squarefree_part(long d);

// For b>0, roundover(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
QUINT roundover(QUINT aa, QUINT bb);

// eclib only has the following for bigints
QUINT sqrt_mod_p(QUINT a, long p); // p odd prime, a quadratic residue

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
