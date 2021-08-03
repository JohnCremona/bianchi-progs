// FILE EUCLID.CC: implementation of pseudo-euclidean algorithm

#define testbezout_psea    // define this to turn on self-verification of bezout via (pseudo-)EA

#include <iostream>

#include "euclid.h"
#include "geometry.h"

/******************************************************************
pseudo-Euclidean step: applies a translation and *if possible* an
M_alpha inversion to a/b (or column vector [a;b]) reducing b, also
multiplying row vector [c.d] my M_alpha on the right.  In the
Euclidean case, the shift is -q where q=a/b (rounded) and the
inversion is via S=[0,-1;1,0], with t=0.  In general if t>=0 then
the t'th inversion was applied.

If the class number is >1 and the ideal (a,b) is non-principal,
then possibly after translation we have that a/b is a singular
point s, in which case no inversion is done and t<0 where s is the
|t|'th singular point.  (The singular points list effectively
starts at index 1.)

a,b,c1,d1,c2,d2 are changed in place, though if either (c2,d2) or
both (c1,d1), (c2,d2) are left as defaults they are not updated.

When applied repeatedly, there are two possible stopping
conditions; note that the ideal (a,b) is unchanged throughout since
we only apply SL(2,O_K)-transformations.  The process is guaranteed
to stop after a finite number of steps since either N(b) is reduced
or the second stopping conditions is reached.

(1) When the ideal (a0,b0) is principal, the stopping condition is
b==0.  Then a = (a0,b0) = a0*d1-b0*d2, and c2/c1=a0/b0 reduced to
lowest terms.

(2) When (a0,b0) is not principal, the stopping condiction is
t<0. Then a/b is the |t|'th singular point, represented as a
fraction with ideal (a,b)=(a0,b0), which may not be the "standard"
representation of the singular point r/s.  Since a/b=r/s, we have
a/r=b/s=lambda, say, where lambda*r=a and lambda*b=s, but *lambda
is not integral* in general.

The quantities c1*a+d1*b andc2*a+d2*b are invariant, i.e. M*v is
invariant where M=[c1,d1;c2,d2] and v=[a;b], since we multuply v on
the left by unimodular matrices at the same time as multiplyin M on
the right by the inverse.

We could use a mat22 [c2,d2;c1,d1] (note the numbering) instead of
c1,d1,c2,d2, but we do not always need one or both pairs: for a gcd
computation, we need neither, for a bezout we need c1,d1 and for
continued fractions we need both.

*********************************************************************/

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
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step with a="<<a<<", N(a)="<<a.nm<<", b="<<b<<", N(b)="<<b.nm<<endl;
#endif
  if (b.nm<0) // impossible unless there has been overflow
    {
      cerr<<"Something is wrong: b="<<b<<" should not have negative norm "<<b.nm<<endl;
      exit(1);
    }
  if (b.nm==0)
    {
      t=0;
      return;
    }

  Quad u, q = 0;  // common simple special case where N(a)<N(b), q = 0 with no work

  if (a.nm<b.nm) // just swap over
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

  q = a/b;  // rounded, so N(a/b -q) is minimal

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

  // Now look for a suitable alpha, trying all in turn (skipping alpha=0)

  Quad r,s, a1,b1;
  mat22 M;
  int local_t = 1;
  for (vector<mat22>::iterator Mi=M_alphas.begin()+1; Mi!=M_alphas.end(); Mi++, local_t++)
    {
      M = *Mi;
      r=-M.d, s=M.c; // alpha = r/s
#ifdef DEBUG_PSEA
      cout<<" - testing type "<<local_t<<", M="<<M<<": ";
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
          cout<<" - success, returning (a,b) = ("<<a<<","<<b<<"), type "<<local_t<<endl;
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

  // We should never arrive here when the class number is 1, as it
  // means that all alphas have failed
  if (Quad::class_number==1)
    {
      cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
      exit(1);
    }

  // Otherwise a/b should be a translate of a singular point sigma, in
  // which case we apply a final translation if necessary and set t=-i
  // where (tranlated) a/b = sigma is the i'th singular point.

  local_t=1;
  RatQuad sigma;
  for (vector<RatQuad>::iterator si=sigmas.begin()+1; si!=sigmas.end(); si++, local_t++)
    {
      sigma = *si;
      r=sigma.num(), s=sigma.den(); // sigma = r/s
#ifdef DEBUG_PSEA
      cout<<" - testing singular point "<<local_t<<": "<<sigma<<flush;
#endif
      // a/b - r/s = (a*s-b*r)/(b*s) will be integral if we have the right sigma
      q = a*s-b*r;
      if (div(b*s,q))  // success!
        {
          q /= (b*s);
          if (q.nm) // else nothing more needs doing since the translation is trivial
            {
              a -= q*b;
              if (compute_c1d1)
                {
                  d1 += q*c1;
                }
              if (compute_c2d2)
                {
                  d2 += q*c2;
                }
            }
#ifdef DEBUG_PSEA
          cout<<" - success, returning (a,b) = ("<<a<<","<<b<<") with a/b = "<<sigma<<" = singular point #"<<local_t<<endl;
#endif
          t = -local_t;
          return;
        }
    } // end of loop over singular points sigma

  // We should never arrive here, as it means that all sigmas have failed
  {
    cerr<<"Pseudo-Euclidean step fails for non-principal ("<<a<<", "<<b<<")"<<endl;
    exit(1);
  }
}

// Declared in quads.h.  Only useful when the ideal (a,b) is
// principal, otherwise 0 is returned.  Can be used to test whether
// (a,b) is principal and return a generator when it is.

Quad quadgcd_psea(const Quad& aa, const Quad& bb)   // Using (pseudo-)EA
{
  if (gcd(aa.nm,bb.nm)==1) return Quad::one;
  Quad a(aa), b(bb); int t=0;
  while (b.nm && t>=0) pseudo_euclidean_step(a, b, t);
  if (b.nm)
    {
      // cout<<"Warning: quadgcd_psea() called with (a,b)=("<<aa<<","<<bb<<"), which generate a non-principal ideal."<<endl;
      // cout<<"  Pseudo-Euclidean Algorithm reached singular point a/b = ("<<a<<")/("<<b<<") = "<<sigmas[-t]<<endl;
      return Quad::zero;
    }
  else
    {
      while (!pos(a)) a*=fundunit;
      return a;
    }
}

// Declared in quads.h.  Only useful when the ideal (a,b) is
// principal, when it returns g such that (g)=(a,b) and x,y such that
// g=x*aa+y*bb. Otherwise 0 is returned, with x, y undefined.

Quad quadbezout_psea(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy)   // Using (pseudo-)EA
{
  Quad a(aa), b(bb), c1(Quad::zero), d1(Quad::one), c2(Quad::one), d2(Quad::zero);
  int t=0;
  while (b.nm>0 && t>=0)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert (c2*d1-c1*d2 == Quad::one);
      assert (c1*a+d1*b==bb);
      assert (c2*a+d2*b==aa);
    }
  if (b.nm)
    {
      // cout<<"Warning: quadbezout_psea() called with (a,b)=("<<aa<<","<<bb<<"), which generate a non-principal ideal."<<endl;
      // cout<<"  Pseudo-Euclidean Algorithm reached singular point a/b = ("<<a<<")/("<<b<<") = "<<sigmas[-t]<<endl;
      return Quad::zero;
    }
  // Now (1) c2*d1-c1*d2 = 1;
  //     (2) c2/c1 = aa/bb (as a reduced fraction), since b = c2*bb-c1*aa = 0;
  //     (3) a = gcd(aa,bb), since (aa,bb)=(a,b)=(a);
  //     (4) aa=a*c2, bb=a*c1;
  //     (5) aa*d1-bb*d2 = a.

  // Note the matrix inversion involved here:
  // [c2,d2;c1,d1]*[a,b] = [aa,bb], so
  // [d1,-d2;-c1,c2]*[aa,bb] = [a,b]

  xx = d1;
  yy = -d2;
  assert (aa*xx+bb*yy==a);

  while (!pos(a))
    {
      a  *= fundunit;
      xx *= fundunit;
      yy *= fundunit;
    }
  if (bb.nm) // reduce x mod bb/g and adjust y to match
    {
      Quad a0 = aa/a, b0 = bb/a;
      Quad t = xx/b0; // rounded
      xx -= b0*t;
      yy += a0*t;
      assert (aa*xx+bb*yy==a);
    }
  return a;
}

// Generalization of extended Euclidean algorithm.

// Given a,b, returns M=[d1,-d2;-c1,c2] (the inverse of [c2,d2;c1,d1])
// such that g/h = M(a/b), i.e. g=d1*a-d2*b, h = -c1*a+c2*b, where
//
// (1) if (a,b) is principal then (a,b)=(g) and h=0, and s=0;
// (2) otherwise (a,b)=(g,h) and g/h is  the s'th singular point (s>=1).
//
// So in Case (1) we get essentially the same information as
// quadbezout_psea(aa, bb, xx, yy) with xx,yy the top row of M.
//
// Note that in case 2, g/h is equal to the singular point sigma=g0/h0
// as an element of the field, but not as a fraction, since the ideal
// (g,h)=(a,b) is in the same (non-principal) ideal class as
// (g0,h0). but is not the same ideal.  In fact, (g,h)=lambda*(g0,h0)
// with lambda = g/g0 = h/h0 (since g*h0=h*g0), but in general lambda
// is not integral.

mat22 generalised_extended_euclid(const Quad& aa, const Quad& bb, int& s)
{
  Quad a(aa), b(bb), c1(Quad::zero), d1(Quad::one), c2(Quad::one), d2(Quad::zero);
  int t=0;
  while (b.nm>0 && t>=0)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert (c2*d1-c1*d2 == Quad::one);
      assert (c1*a+d1*b==bb);
      assert (c2*a+d2*b==aa);
    }

  // Now (1) c2*d1-c1*d2 = 1;
  //     (2) c2/c1 = aa/bb (as a reduced fraction), since b = c2*bb-c1*aa = 0;
  //     (3) a = gcd(aa,bb), since (aa,bb)=(a,b)=(a);
  //     (4) aa=a*c2, bb=a*c1;
  //     (5) aa*d1-bb*d2 = a.

  // Note the matrix inversion involved here:
  // [c2,d2;c1,d1]*[a,b] = [aa,bb], so
  // [d1,-d2;-c1,c2]*[aa,bb] = [a,b]

  s = (b.nm==0? 0 : -t);

  mat22 M(d1,-d2,-c1,c2); // maps aa/bb to a/b
  assert (d1*aa-d2*bb == a);
  assert (-c1*aa+c2*bb == b);
  return M;
}
