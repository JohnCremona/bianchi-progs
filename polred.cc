// FILE POLRED.CC: implementation of functions for reducing ZZX polynomials (polredabs) via libpari

#include "polred.h"
#include <eclib/interface.h> // for getenv_with_default
#include <eclib/pari_init.h>
#include <eclib/convert.h>
#include <assert.h>

using PARI::degpol;
using PARI::t_POL;
using PARI::RgX_renormalize_lg;
using PARI::polredabs;
using PARI::polredabs0;
using PARI::nf_ORIG;
using PARI::lift;

//#define DEBUG_POLY

// convert a t_POL (with Z coefficients) to a ZZX
ZZX t_POL_to_ZZX(GEN P)
{
#ifdef DEBUG_POLY
  pari_printf(" Converting t_POL %Ps to ZZX\n", P);
#endif
  int d = degpol(P);
  ZZX f;
  for (int i=0; i<=d; i++)
    SetCoeff(f, i, PARI_to_NTL(gel(P, i+2)));
#ifdef DEBUG_POLY
  cout << " Result is " << str(f) << endl;
#endif
  return f;
}

GEN ZZX_to_t_POL(const ZZX& f)
{
#ifdef DEBUG_POLY
  cout << " Converting  ZZX " << str(f) << " to t_POL" << endl;
#endif
  int d = deg(f);
  GEN P = cgetg(d+3, t_POL);
  P[1] = evalvarn(0); // set variable to #0, i.e. 'x'
  for (int i=0; i<=d; i++)
    gel(P,i+2) = NTL_to_PARI(coeff(f,i));
  P = RgX_renormalize_lg(P,d+3);
#ifdef DEBUG_POLY
  pari_printf(" Result is %Ps\n", P);
#endif
  return P;
}

// polredabs of an *irreducible* polynomial in Z[X]
// (1) return monic integral g defining the same field as f
ZZX polredabs(const ZZX& f)
{
  if (!IsIrreducible(f))
    {
      cerr << "polredabs() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN G = polredabs(ZZX_to_t_POL(f));
  pari_printf("polredabs(f) returns %Ps\n", G);
  ZZX g = t_POL_to_ZZX(G);
  avma = av;
  return g;
}

// (2) also sets h such that a=h(b) (so f(h(b))=0) where f(a)=g(b)=0

ZZX polredabs(const ZZX& f, ZZX& h)
{
  if (!IsIrreducible(f))
    {
      cerr << "polredabs() called with f = " << str(f) << " which is reducible" << endl;
      return f;
    }
  pari_sp av = avma;
  GEN G_H = polredabs0(ZZX_to_t_POL(f), nf_ORIG);
  pari_printf("polredabs0(f,nf_ORIG) returns [G, H] = %Ps\n", G_H);
  GEN G = gel(G_H,1);
  pari_printf("  G = %Ps\n", G);
  GEN H = lift(gel(G_H,2));
  pari_printf("  H = %Ps\n", H);
  ZZX g = t_POL_to_ZZX(G);
  h = t_POL_to_ZZX(H);
  avma = av;
  return g;
}

