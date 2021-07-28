// FILE EUCLID.H: declaration of pseudo-euclidean algorithm function

#if     !defined(_EUCLID_H)
#define _EUCLID_H      1       //flags that this file has been included

#include <iostream>

#include "quads.h"

// pseudo-Euclidean step: applies a translation and *if possible* an
// M_alpha inversion to a/b (or column vector [a;b]) reducing b, also
// multiplying row vector [c.d] my M_alpha on the right.  In the
// Euclidean case, the shift is -q where q=a/b (rounded) and the
// inversion is via S=[0,-1;1,0], with t=0.  In general if t>=0 then
// the t'th inversion was applied.

// If the class number is >1 and the ideal (a,b) is non-principal,
// then possibly after translation we have that a/b is a singular
// point s, in which case no inversion is done and t<0 where s is the
// |t|'th singular point.  (The singular points list effectively
// starts at index 1.)

// a,b,c1,d1,c2,d2 are changed in place, though if either (c2,d2) or
// both (c1,d1), (c2,d2) are left as defaults they are not updated.

// When applied repeatedly, there are two possible stopping
// conditions; note that the ideal (a,b) is unchanged throughout since
// we only apply SL(2,O_K)-transformations.  The process is guaranteed
// to stop after a finite number of steps since either N(b) is reduced
// or the second stopping conditions is reached.

// (1) When the ideal (a0,b0) is principal, the stopping condition is
// b==0.  Then a = (a0,b0) = a0*d1-b0*d2, and c2/c1=a0/b0 reduced to
// lowest terms.

// (2) When (a0,b0) is not principal, the stopping condiction is
// t<0. Then a/b is the |t|'th singular point, represented as a
// fraction with ideal (a,b)=(a0,b0), which may not be the "standard"
// representation of the singular point r/s.  Since a/b=r/s, we have
// a/r=b/s=lambda, say, where lambda*r=a and lambda*b=s, but *lambda
// is not integral* in general.


void pseudo_euclidean_step(Quad& a, Quad& b, int& t,
                           Quad& c1=Quad::zero, Quad& d1=Quad::zero,
                           Quad& c2=Quad::zero, Quad& d2=Quad::zero);

#endif
