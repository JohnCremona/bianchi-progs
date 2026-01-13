// FILE POLRED.H: declaration of functions for reducing ZZX polynomials (polredabs) via libpari

#if     !defined(_POLRED_H)
#define _POLRED_H      1       //flags that this file has been included

#include <eclib/templates.h>
#include <eclib/int.h>
#include <eclib/polys.h>
#include <eclib/pari_init.h>
#undef recip // pariold.h #defines recip = serreverse

// conversion from ZZX to t_POL and back:
ZZX t_POL_to_ZZX(GEN P);
GEN ZZX_to_t_POL(const ZZX& f);

// polredabs of an *irreducible* polynomial in Z[X]
// (1) return monic integral g defining the same field as f
ZZX polredabs(const ZZX& f);
// (2) also sets h such that a=h(b) (so f(h(b))=0) where f(a)=g(b)=0
ZZX polredabs(const ZZX& f, ZZX& h);
#endif
