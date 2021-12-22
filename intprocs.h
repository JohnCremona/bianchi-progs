// FILE INTPROCS.H: miscellaneous integer (int/long) functions

#if     !defined(_INTPROCS_H)
#define _INTPROCS_H      1       //flags that this file has been included

#include <eclib/interface.h>
#include <eclib/marith.h>

typedef bigint QUINT; // integer type for components of a Quad

int is_nonnegative(QUINT a);
long squarefree_part(long d);

// For b>0, rounded_division(a,b) = q such that a/b = q + r/b with -1/2 <= r/b < 1/2
//
QUINT rounded_division(QUINT aa, QUINT bb);
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
