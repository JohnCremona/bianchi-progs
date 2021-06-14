// FILE EUCLID.CC: implementation of pseudo-euclidean algorithm

#define testbezout_psea    // define this to turn on self-verification of bezout via (pseudo-)EA

#include <iostream>

#include <eclib/arith.h>
#include "quads.h"
#include "euclid.h"
#include "geometry.h"

// pseudo-Euclidean step: applies a translation and M_alpha
// to a/b (or column vector [a;b]) reducing b, also multiplying row
// vectors [c1,d1] and [c2,d2] by M_alpha^{-1} on the right.  In the Euclidean case, the
// shift is -q where q=a/b (rounded) and the inversion is via S=[0,-1;1,0] (alpha=0).

// a,b,c,d are changed in place, and on return, t holds the "type"
// (index of alpha which worked)

//#define DEBUG_PSEA

void pseudo_euclidean_step(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2)
{
  // We update c1,d1 unless they are both 0 and similarly c2,d2.
  // We record the type of the transformation in t unless it is initialised to -1.
  // For simple gcd, we need none of these;  for bezout (extended EA) we need c1,d1 and c2,d2
  // For convergents we need c1,d1 and t

  int compute_c1d1=1, compute_c2d2=1;
  if (c1==Quad::zero && d1==Quad::zero)
    compute_c1d1 = 0;
  if (c2==Quad::zero && d2==Quad::zero)
    compute_c2d2 = 0;

  if (b.nm==0)
    {
      t=0;
      return;
    }
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step("<<a<<","<<b<<"), N(b)="<<b.nm<<endl;
#endif

  Quad u, q = 0;  // common simple special case where N(a)<N(b), q = 0 with no work
  if (a.nm >= b.nm)
    {
      q = a/b;  // rounded, so N(a/b -q) is minimal
    }

#ifdef DEBUG_PSEA
  cout<<" - translation = "<<q<<endl;
#endif
  a -= q*b;
  if (compute_c1d1) d1 += q*c1;
  if (compute_c2d2) d2 += q*c2;
#ifdef DEBUG_PSEA
  cout<<" - reduced a = "<<a<<endl;
#endif
  if (a.nm < b.nm) // always true in Euclidean case; invert using S (N.B. S^{-1}=-S)
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0";
      cout << endl;
#endif
      t = 0;
      return;
    }

  // Now look or a suitable alpha, trying all in turn (skipping alpha=0)

  Quad r,s, a1,b1;
  mat22 M;
  int local_t = 1;
  for (vector<mat22>::iterator Mi=M_alphas.begin()+1; Mi!=M_alphas.end(); Mi++, local_t++)
    {
      M = *Mi;
      r=-M.d, s=M.c; // alpha = r/s
#ifdef DEBUG_PSEA
      cout<<" - testing type "<<local_t<<", M="<<M<<endl;
#endif

      // We need to use temporary copies of a,b in case this alpha fails
      a1 = a; b1 = b; q = 0;

      // First see if we can reduce without a further shift (since
      // rounded division is relatively expensive, and here the new
      // quotient q will usually be 0 anyway).

      M.apply_left(a1,b1);
      if (b1.nm >= b.nm) // not successful yet
        {
          // Find the shift taking a/b closest to alpha
          q = (a*s-b*r)/(b*s); // closest integer to (a/b)-(r/s)
          if (q.nm)            // do the extra translation by q
            {
              a1 = a-q*b, b1 = b;
              M.apply_left(a1,b1);
            }
        }
      if (b1.nm < b.nm) // success!
        {
          a = a1;
          b = b1;
          if (compute_c1d1)
            {
              if (q.nm) d1 += q*c1;
              M.apply_right_inverse(c1,d1);
            }
          if (compute_c2d2)
            {
              if (q.nm) d2 += q*c2;
              M.apply_right_inverse(c2,d2);
            }
#ifdef DEBUG_PSEA
          cout<<" - success, returning (a,b) = ("<<a<<","<<b<<")";
          cout <<", type "<<local_t;
          cout<<endl;
#endif
          t = local_t;
          return;
        }
#ifdef DEBUG_PSEA
      else
        {
          cout<<" - failure (q="<<q<<"), new b would have had norm "<<quadnorm(b1)<<endl;
        }
#endif
    }
  // We should never arrive here as it means that all alphas have failed
  t = -1;
  cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
  exit(1);
}

Quad quadgcd_psea(const Quad& aa, const Quad& bb)   // Using (pseudo-)EA
{
  if (gcd(aa.nm,bb.nm)==1) return Quad::one;
  Quad a(aa), b(bb); int t;
  while (b.nm) pseudo_euclidean_step(a, b, t);
  while (!pos(a)) a*=fundunit;
  return a;
}

Quad quadbezout_psea(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy)   // Using (pseudo-)EA
{
  long x,y;
  if (bezout(aa.nm,bb.nm,x,y)==1)
    {
      xx=x*quadconj(aa);
      yy=y*quadconj(bb);
      if (bb.nm)
        {
          Quad t = xx/bb; // rounded
          xx -= bb*t;
          yy += aa*t;
        }
      assert (aa*xx+bb*yy==Quad::one);
      return Quad::one;
    }
  Quad a(aa), b(bb), c1(Quad::zero), d1(Quad::one), c2(Quad::one), d2(Quad::zero);
  int t;
  while (b.nm)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert (c2*d1-c1*d2 == Quad::one);
      assert (c1*a+d1*b==bb);
      assert (c2*a+d2*b==aa);
      assert (c2*d1-c1*d2 == Quad::one);
    }
  // Now (1) c2*d1-c1*d2 = 1;
  //     (2) c2/c1 = aa/bb (as a reduced fraction);
  //     (3) a = gcd(aa,bb);
  //     (4) aa=a*c2, bb=a*c1;
  //     (5) aa*d1-bb*d2 = a.
  assert (aa*d1-bb*d2 == a);
  xx = d1;
  yy = -d2;
  while (!pos(a))
    {
      a  *= fundunit;
      xx *= fundunit;
      yy *= fundunit;
    }
  if (bb.nm)
    {
      Quad a0 = aa/a, b0 = bb/a;
      Quad t = xx/b0; // rounded
      xx -= b0*t;
      yy += a0*t;
    }
  assert (aa*xx+bb*yy==a);
#ifdef testbezout_psea
//CHECK:
  if (div(a,aa) && div(a,bb) && (a==xx*aa+yy*bb)) {;}  //OK
  else
    {cerr<<"Error in quadbezout_psea!"<<endl;
     cerr<<"a = "<<aa<<endl;
     cerr<<"b = "<<bb<<endl;
     cerr<<"x = "<<xx<<endl;
     cerr<<"y = "<<yy<<endl;
     cerr<<"g = "<<a<<endl;
    }
#endif
  return a;
}
