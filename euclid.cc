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

// This function is crucial in reducing modular symbols, expressing
// each as a linear combination of generalised M-symbols.

//#define DEBUG_PSEA

void pseudo_euclidean_step_orig(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2);
void pseudo_euclidean_step_old(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2);
void pseudo_euclidean_step_new(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2);
void pseudo_euclidean_step(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2)
{
  pseudo_euclidean_step_old(a, b, t, c1, d1, c2,  d2);
}

//return type t if translation by -q followed by M_alphas[t] reduces
//(a,b) to (a',b') with N(b')<N(b), or -1 if none exists, which will
//be if and only if a/b is singular.  If lucky, return the first t
//which reduces N(b), otherwise try all and return the best.
int find_best_alpha_old(const Quad& a, const Quad& b, Quad& shift, int lucky=0);
int find_best_alpha_new(const Quad& a, const Quad& b, Quad& shift, int lucky=0);
int find_best_alpha(const Quad& a, const Quad& b, Quad& shift, int lucky=1)
{
  return find_best_alpha_old(a, b, shift, lucky);
}

// This assumes that a is reduced mod b and that N(a)>N(b) so alpha=0 will not work
int find_best_alpha_old(const Quad& a, const Quad& b, Quad& shift, int lucky)
{
  Quad a0 = a, b0=b, b1, s, best_shift(0);
  INT best_bnorm = b.norm(), b0norm=b0.norm();
  int t=0, best_t=-1;
  // Look for a suitable alpha, trying all in turn, returning the type
  // of the one which gives best reduction
  for (const auto& M: M_alphas)
    {
      if (t==0) // skipping alpha=0
        {
          t++;
          continue;
        }
      // First see if we can reduce without a further shift
      shift = 0;
      b1 = M.apply_left_den(a0,b0);
      if (b1.norm() >= b0norm) // not successful yet
        {
          // Find the shift taking a/b closest to alpha
          s = M.entry(1,0);
          shift = b1/(b0*s);      // closest integer to (a1/b1)-(r/s)
          if (!shift.is_zero())  // do the extra translation by q
            b1 = M.apply_left_den(a0-shift*b0,b0);
        }
      if (b1.norm() < best_bnorm)
        {
          best_bnorm = b1.norm();
          best_shift = shift;
          best_t     = t;
          if (lucky)
            break;
        }
      t++;
    }
  shift = best_shift;
  return best_t;
}

int find_best_alpha_new(const Quad& a, const Quad& b, Quad& shift, int lucky)
{
  // (1) Look for a suitable alpha = r/s: loop over s in alpha_denoms
  // and see if there exists r such that N((a/b)*s-r)<1, with (r,s)
  // principal.  If so, r/s will be a translate of an alpha which
  // works.  Keep going until the best one is found.
  // First reduce a/b to rectangle
  Quad best_shift(0);
  INT best_bnorm = b.norm();
  int best_t=-1;
  for ( const auto& s : alpha_denoms)
    {
      //auto rlist1 = nearest_quads(RatQuad(s*a, b), 0); // 0 means all (there may be two), 1 means at most one
      Quad as = a*s;
      auto rlist = nearest_quads_to_quotient(as, b, 0); // 0 means all (there may be two), 1 means at most one
      for (const auto& r : rlist) // loop does nothing of rlist is empty
        {
          int t = alpha_index_with_translation(RatQuad(r,s), shift);
          if (t==-1)
            continue;
          // Apply the translation and check that M does now reduce (a,b)
          Quad b1 = M_alphas[t].apply_left_den(a-shift*b,b);
          if (b1.norm() < best_bnorm)
            {
              best_bnorm = b1.norm();
              best_shift = shift;
              best_t     = t;
              if (lucky)
                break;
            }
        } // end of r loop
    } // end of s loop
  shift = best_shift;
  return best_t;
}

//#undef DEBUG_PSEA

void pseudo_euclidean_step_new(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2)
{
  // We update c1,d1 unless they are both 0 and similarly c2,d2.
  // We record the type of the transformation in t unless it is initialised to -1.
  // For simple gcd, we need none of these;  for bezout (extended EA) we need c1,d1 and c2,d2
  // For convergents we need c1,d1 and t
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step with a="<<a<<", N(a)="<<a.norm()<<", b="<<b<<", N(b)="<<b.norm()<<endl;
#endif
  t = 0;
  if (b.is_zero())
    return;

  int compute_c1d1 = !(c1.is_zero() && d1.is_zero());
  int compute_c2d2 = !(c2.is_zero() && d2.is_zero());
  Quad original_a = a, original_b = b, q, u;

  // (0) common easy special case where N(a)<N(b) after a translation (type 0)

  q = a/b; // rounded
#ifdef DEBUG_PSEA
  cout<<" - initial translation = "<<q<<endl;
#endif
  a.subprod(q,b);
  if (compute_c1d1) d1.addprod(q,c1);
  if (compute_c2d2) d2.addprod(q,c2);

  if (a.norm() < b.norm()) // always true in Euclidean case; invert using S (N.B. S^{-1}=-S)
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0" << endl;
#endif
      return;
    }

  // Now the real work begins.

  // (1) Look for a suitable alpha

  t = find_best_alpha_new(a, b, q, 1); // 1 means use the first one found
  if (t!=-1) // none found
    {
      mat22 M = M_alphas[t];
      a.subprod(q,b);
      M.apply_left(a,b);
      if (compute_c1d1)
        {
          d1.addprod(q,c1);
          M.apply_right_inverse(c1,d1);
        }
      if (compute_c2d2)
        {
          d2.addprod(q,c2);
          M.apply_right_inverse(c2,d2);
        }
#ifdef DEBUG_PSEA
      cout<<" - reduction success (with shift="<<q<<" and alpha #"<<t
          <<"), returning (a,b) = ("<<a<<","<<b<<"), type "<<t<<endl;
#endif
      return;
    }

  // (2) Otherwise, a/b is a translate of a singular point sigma, we apply an
  // extra translation if necessary and set t=-i where the translate
  // of a/b is the i'th singular point.  This is quicker to check so
  // we do it first; if this test fails, then a/b can be reduced using
  // at least one alpha.

#ifdef DEBUG_PSEA
  cout<<" - no reducing alpha found, a/b must be singular"<<endl;
#endif

  t = sigma_index_with_translation(a, b, q);
  if (t==-1)
    {
      // We should never arrive here, as it means that all alphas have
      // failed though a/b is not singular.
      cerr<<"Pseudo-Euclidean step fails for ("<<original_a<<", "<<original_b<<")"<<" in field d="<<Quad::d<<endl;
      exit(1);
    }
  a.subprod(q,b);
  if (compute_c1d1) d1.addprod(q,c1);
  if (compute_c2d2) d2.addprod(q,c2);
#ifdef DEBUG_PSEA
  cout<<" - success, returning (a,b) = ("<<a<<","<<b<<") with a/b = " << sigmas[t]
      <<" = singular point #"<<t<<endl;
#endif
  t = -t;
  return;
}

void pseudo_euclidean_step_old(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2)
{
  // We update c1,d1 unless they are both 0 and similarly c2,d2.
  // We record the type of the transformation in t unless it is initialised to -1.
  // For simple gcd, we need none of these;  for bezout (extended EA) we need c1,d1 and c2,d2
  // For convergents we need c1,d1 and t
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step with a="<<a<<", N(a)="<<a.norm()<<", b="<<b<<", N(b)="<<b.norm()<<endl;
#endif
  t = 0;
  if (b.is_zero())
    return;

  int compute_c1d1 = !(c1.is_zero() && d1.is_zero());
  int compute_c2d2 = !(c2.is_zero() && d2.is_zero());
  Quad u, shift = a/b;  // rounded, so N(a/b - shift) is minimal

  a.subprod(shift,b);
  if (compute_c1d1) d1.addprod(shift,c1);
  if (compute_c2d2) d2.addprod(shift,c2);
#ifdef DEBUG_PSEA
  cout<<" - translation = "<<shift<<endl;
#endif

  if (a.norm()<b.norm()) // standard Euclidean M works
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0" << endl;
#endif
      t = 0;
      return;
    }

  // If a/b is a translate of a singular point sigma, we apply an
  // extra translation if necessary and set t=-i where the translate
  // of a/b is the i'th singular point.  This is quick to check so we
  // do it first; if this fails then a/b can be reduced using at least
  // one alpha.

  t = sigma_index_with_translation(a,b, shift);
  if (t!=-1)
    {
      a.subprod(shift,b);
      if (compute_c1d1) d1.addprod(shift,c1);
      if (compute_c2d2) d2.addprod(shift,c2);
#ifdef DEBUG_PSEA
      cout<<" - a/b is singular, a translate of sigma["<<t<<"]"<<endl;
#endif
      t = -t;
      return;
    }

  // Now look for a suitable alpha (skipping 0/1 which we already tested)

  t = find_best_alpha_old(a, b, shift, 1); // 1 means use the first one found, 0 find the best
  if (t==-1)
    {
#ifdef DEBUG_PSEA
      cout<<" - all alphas failed thoough a/b is not singular"<<endl;
#endif
      // We should never arrive here, as it means that all alphas have
      // failed and a/b is not singular.
      cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
      exit(1);
    }

  mat22 M = M_alphas[t];
  a.subprod(shift,b);
  M.apply_left(a,b);
  if (compute_c1d1)
    {
      d1.addprod(shift,c1);
      M.apply_right_inverse(c1,d1);
    }
  if (compute_c2d2)
    {
      d2.addprod(shift,c2);
      M.apply_right_inverse(c2,d2);
    }
#ifdef DEBUG_PSEA
  cout<<" - reduction success (with shift="<<shift<<" and alpha="<<alphas[t]
      <<"), returning (a,b) = ("<<a<<","<<b<<"), type "<<t<<endl;
#endif
  return;
}

void pseudo_euclidean_step_orig(Quad& a, Quad& b, int& t, Quad& c1, Quad& d1, Quad& c2, Quad& d2)
{
  // We update c1,d1 unless they are both 0 and similarly c2,d2.
  // We record the type of the transformation in t unless it is initialised to -1.
  // For simple gcd, we need none of these;  for bezout (extended EA) we need c1,d1 and c2,d2
  // For convergents we need c1,d1 and t
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step with a="<<a<<", N(a)="<<a.norm()<<", b="<<b<<", N(b)="<<b.norm()<<endl;
#endif
  if (b.is_zero())
    {
      t=0;
      return;
    }

  int compute_c1d1 = !(c1.is_zero() && d1.is_zero());
  int compute_c2d2 = !(c2.is_zero() && d2.is_zero());
  Quad u, q = Quad::zero;  // common simple special case where N(a)<N(b), q = 0 with no work

  if (a.norm()<b.norm()) // just swap over
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0" << endl;
#endif
      t = 0;
      return;
    }

  q = a/b;  // rounded, so N(a/b -q) is minimal

#ifdef DEBUG_PSEA
  cout<<" - translation = "<<q<<endl;
#endif
  if (!q.is_zero())
    {
      a -= q*b;
      if (compute_c1d1) d1 += q*c1;
      if (compute_c2d2) d2 += q*c2;
    }
#ifdef DEBUG_PSEA
  cout<<" - reduced a = "<<a<<endl;
#endif
  if (a.norm() < b.norm()) // always true in Euclidean case; invert using S (N.B. S^{-1}=-S)
    {
      u = a; a=-b; b=u;
      if (compute_c1d1) {u = -d1; d1=c1; c1=u;}
      if (compute_c2d2) {u = -d2; d2=c2; c2=u;}
#ifdef DEBUG_PSEA
      cout<<" - after inverting by S, returning (a,b) = ("<<a<<","<<b<<") ";
      if (compute_c1d1) cout << "(c1,d1)=("<<c1<<","<<d1<<") ";
      if (compute_c2d2) cout << "(c2,d2)=("<<c2<<","<<d2<<") ";
      cout <<" type=0" << endl;
#endif
      t = 0;
      return;
    }

  // If a/b is a translate of a singular point sigma, we apply an
  // extra translation if necessary and set t=-i where the translate
  // of a/b is the i'th singular point.  This is quick to check so we
  // do it first; if this fails then a/b can be reduced using at least
  // one alpha.

  Quad r,s,bs;
  int local_t=1;
  for (const auto& sig : sigmas)
    {
      if (sig.is_infinity()) continue;
#ifdef DEBUG_PSEA
      cout<<" - testing singular point "<<local_t<<": "<< sig <<endl;
#endif
      r=sig.num(), s=sig.den(); // sigma = r/s
      // a/b - r/s = (a*s-b*r)/(b*s) will be integral if we have the right sigma
      q = a*s-b*r;
      bs = b*s;
      if (div(bs,q))  // success!
        {
          q /= bs; // rounded
          if (!q.is_zero()) // else nothing more needs doing since the translation is trivial
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
          cout<<" - success, returning (a,b) = ("<<a<<","<<b<<") with a/b = "<< sig
              <<" = singular point #"<<local_t<<endl;
#endif
          t = -local_t;
          return;
        }
      local_t++;
    } // end of loop over singular points sigma

#ifdef DEBUG_PSEA
  cout<<" - a/b is not singular, looking for an alpha which reduces"<<endl;
#endif

  // Now look for a suitable alpha, trying all in turn (skipping alpha=0 which we checked already)

  Quad a1,b1;
  local_t = 1;
  int first = 1;
  for ( auto& M : M_alphas)
    {
      if (first) {first=0;continue;} // so we skip M_alphas[0]
      r=-M.entry(1,1), s=M.entry(1,0); // alpha = r/s
#ifdef DEBUG_PSEA
      //cout<<" - testing type "<<local_t<<", M="<<M<<": ";
#endif

      // We need to use temporary copies of a,b in case this alpha fails
      a1 = a; b1 = b; q = Quad::zero;

      // First see if we can reduce without a further shift (since
      // rounded division is relatively expensive, and here the new
      // quotient q will usually be 0 anyway).

      M.apply_left(a1,b1);
      if (b1.norm() >= b.norm()) // not successful yet
        {
          // Find the shift taking a/b closest to alpha
          q = (a*s-b*r)/(b*s); // closest integer to (a/b)-(r/s)
          if (!q.is_zero())            // do the extra translation by q
            {
              a1 = a-q*b, b1 = b;
              M.apply_left(a1,b1);
            }
        }
      if (b1.norm() < b.norm()) // success!
        {
          a = a1;
          b = b1;
          if (compute_c1d1)
            {
              if (!q.is_zero()) d1 += q*c1;
              M.apply_right_inverse(c1,d1);
            }
          if (compute_c2d2)
            {
              if (!q.is_zero()) d2 += q*c2;
              M.apply_right_inverse(c2,d2);
            }
#ifdef DEBUG_PSEA
          cout<<" - success (q="<<q<<"), returning (a,b) = ("<<a<<","<<b<<"), type "<<local_t<<endl;
#endif
          t = local_t;
          return;
        }
#ifdef DEBUG_PSEA
      else
        {
          //cout<<" - failure (q="<<q<<"), new b would have had norm "<<b1.norm()<<endl;
        }
#endif
      local_t++;
    }

  // We should never arrive here, as it means that all alphas have
  // failed and a/b is not singular.
  cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
  exit(1);
}


// The next two functions (declared in quads.h) are not used, as the
// default versions quadgcd_default and quadbezout_default using
// ideals are simpler and faster.

// Declared in quads.h.  Only useful when the ideal (a,b) is
// principal, otherwise 0 is returned.  Can be used to test whether
// (a,b) is principal and return a generator when it is.

// NB using quadgcd_default (using ideals) is simpler and probably faster

Quad quadgcd_psea(const Quad& aa, const Quad& bb)   // Using (pseudo-)EA
{
  if (gcd(aa.norm(),bb.norm()).is_one()) return Quad::one;
  Quad a(aa), b(bb); int t=0;
  while (!b.is_zero() && t>=0) pseudo_euclidean_step(a, b, t);
  if (!b.is_zero())
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
  // cout<<"quadbezout("<<aa<<","<<bb<<")"<<endl;
  Quad a(aa), b(bb), c1(Quad::zero), d1(Quad::one), c2(Quad::one), d2(Quad::zero);
  int t=0;
  while (b.norm()>0 && t>=0)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert ((c2*d1-c1*d2).is_one());
      assert (c1*a+d1*b==bb);
      assert (c2*a+d2*b==aa);
    }
  if (!b.is_zero())
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
  // cout<<"quadbezout("<<aa<<","<<bb<<") finds gcd="<<a<<" and initially xx="<<xx<<", yy="<<yy<<endl;

  while (!pos(a))
    {
      a  *= fundunit;
      xx *= fundunit;
      yy *= fundunit;
    }
  if (!bb.is_zero()) // reduce x mod bb/g and adjust y to match
    {
      Quad a0 = aa/a, b0 = bb/a;
      Quad tmp = xx/b0; // rounded
      // cout<<"quadbezout("<<aa<<","<<bb<<") next finds gcd="<<a<<" and xx="<<xx<<", yy="<<yy<<endl;
      // cout<<" adjustment t="<<tmp<<" is rounded division of xx by "<<bb<<"/"<<a<<endl;
      xx -= b0*tmp;
      yy += a0*tmp;
      // cout<<" finally xx="<<xx<<", yy="<<yy<<endl;
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
  while (b.norm()>0 && t>=0)
    {
      pseudo_euclidean_step(a, b, t, c1, d1, c2, d2);
      assert ((c2*d1-c1*d2).is_one());
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

  s = (b.norm().is_zero()? 0 : -t);

  mat22 M(d1,-d2,-c1,c2); // maps aa/bb to a/b
  assert (d1*aa-d2*bb == a);
  assert (-c1*aa+c2*bb == b);
  return M;
}
