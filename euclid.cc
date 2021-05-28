// FILE EUCLID.CC: implementation of pseudo-euclidean algorithm and associated data

#define testbezout_psea    // define this to turn on self-verification of bezout via (pseudo-)EA

#include <iostream>

#include <eclib/arith.h>
#include "quads.h"
#include "euclid.h"

// Definitions of commonly used matrices

mat22 mat22::identity(1,0,0,1);
mat22 mat22::J(-1,0,0,1);
mat22 mat22::S(0,-1,1,0);
mat22 mat22::T(1,1,0,1);
mat22 mat22::U(1,Quad::w,0,1);
mat22 mat22::TS(1,-1,1,0);   // = T*S
mat22 mat22::TiS(-1,-1,1,0); // = T^{-1}*S
mat22 mat22::R(0,1,1,0);

// Definitions of alphas and associated matrices M_alpha such that
// det(M_alpha)=1 and M_alpha(alpha)=oo.
//
// alpha_pairs is a permutation of range(N_alphas) such that
// alpha_pairs[i]=j where M_alpha[i](oo) = alpha[j].
//
// We do not store the alphas explicitly, they are -d/c where M_alpha=[a,b;c,d].

int n_alphas;
vector<mat22> M_alphas;
vector<int> alpha_pairs;

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  mat22 M_alpha(a,b,c,d);  // maps alpha = -d/c to oo
  assert (M_alpha.det()==1);
  M_alphas.push_back(M_alpha);
  n_alphas += 1;
}

// Global function to be used once after setting the field:

void define_alphas()
{
  int d = Quad::d;

  // alphas (only 0) with denominator 1:

  add_alpha(0,-1,1,0);  // alpha[0] = 0
  alpha_pairs.push_back(0); // 0-0

  if (d<19) return;

  Quad w = Quad::w;

  // alphas (w/2, (w-1)/2) with denominator 2:

  Quad u = (d-3)/8;  // = 2, 5, 8, 20
  add_alpha(w-1,u,2,-w);  // alpha[1] = w/2
  add_alpha(w,u,2,1-w);   // alpha[2] = (w-1)/2
  alpha_pairs.push_back(2); // 1-2
  alpha_pairs.push_back(1); // 2-1

  if (d<43) return;

  // alphas (w/3, -w/3, (1-w)/3, (w-1)/3, (w+1)/3, -(w+1)/3) with denominator 3:

  u = -(d+5)/12;  // = -4, -6, -14 so w^2 = w+3*u+1
  add_alpha(1-w,u,3,-w);           // alpha[3] = w/3
  add_alpha(w-1,u,3,w);            // alpha[4] = -w/3
  add_alpha(w,u,3,w-1);            // alpha[5] = (1-w)/3
  add_alpha(-w,u,3,1-w);           // alpha[6] = (w-1)/3
  add_alpha(w+1,-(w+u+1),3,-1-w);  // alpha[7] = (w+1)/3
  add_alpha(-w-1,-(w+u+1),3,1+w);  // alpha[8] = -(w+1)/3
  alpha_pairs.push_back(5); // 3-5
  alpha_pairs.push_back(6); // 4-6
  alpha_pairs.push_back(3); // 5-3
  alpha_pairs.push_back(4); // 6-4
  alpha_pairs.push_back(7); // 7-7
  alpha_pairs.push_back(8); // 8-8

  if (d<67) return;

  // alphas with denominator 4 for both d=67 and d=163:

  u = (d-3)/16; // = 4, 10 so w^2 = w - (4*u+1)
  add_alpha(w-1, u, 4, -w); // alpha[9] = w/4
  add_alpha(1-w, u, 4, w);  // alpha[10] = -w/4
  add_alpha(w, u, 4, 1-w);  // alpha[11] = (w-1)/4
  add_alpha(-w, u, 4, w-1); // alpha[12] = (1-w)/4
  alpha_pairs.push_back(11); // 9-11
  alpha_pairs.push_back(12); // 10-12
  alpha_pairs.push_back(9);  // 11-9
  alpha_pairs.push_back(10); // 12-10

  add_alpha(2-w, -u-1, 4, -w-1); // alpha[13] = (1+w)/4
  add_alpha(w-2, -u-1, 4, w+1);  // alpha[14] = -(1+w)/4
  add_alpha(w+1, -u-1, 4, w-2);  // alpha[15] = (2-w)/4
  add_alpha(-1-w,-u-1, 4, 2-w);  // alpha[16] = (w-2)/4
  alpha_pairs.push_back(15);     // 13-15
  alpha_pairs.push_back(16);     // 14-16
  alpha_pairs.push_back(13);     // 15-13
  alpha_pairs.push_back(14);     // 16-14

  if (d==67)   // alphas with denominator norm 23 for d=67 only:
    {
      Quad s = 3-w; // norm 23
      add_alpha(2+w, 7-w, s, -w-6); // alpha[17] = (6+w)/(3-w)
      add_alpha(-2-w, 7-w, s, w+6); // alpha[18] = (-6-w)/(3-w)
      add_alpha(6+w, 7-w, s, -w-2); // alpha[19] = (2+w)/(3-w)
      add_alpha(-6-w, 7-w, s, w+2); // alpha[20] = (-2-w)/(3-w)
      alpha_pairs.push_back(19);     // 17-19
      alpha_pairs.push_back(20);     // 18-20
      alpha_pairs.push_back(17);     // 19-17
      alpha_pairs.push_back(18);     // 20-18

      s = 2+w; // norm 23, conjugate to previous
      add_alpha(3-w, 6+w, s, w-7); // alpha[21] = (7-w)/(2+w)
      add_alpha(w-3, 6+w, s, 7-w); // alpha[22] = (w-7)/(2+w)
      add_alpha(7-w, 6+w, s, w-3); // alpha[23] = (3-w)/(2+w)
      add_alpha(w-7, 6+w, s, 3-w); // alpha[24] = (w-3)/(2+w)
      alpha_pairs.push_back(23);     // 21-23
      alpha_pairs.push_back(24);     // 22-24
      alpha_pairs.push_back(21);     // 23-21
      alpha_pairs.push_back(22);     // 24-22

      return;
    }
  cerr << "define_alphas() not yet implemented for field "<<d<<endl;
}


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

  int compute_c1d1=1, compute_c2d2=1, compute_t=1;
  if (c1==Quad::zero && d1==Quad::zero)
    compute_c1d1 = 0;
  if (c2==Quad::zero && d2==Quad::zero)
    compute_c2d2 = 0;
  if (t==default_t)
    compute_t = 0;

  if (b.nm==0)
    {
      if (compute_t) t=0;
      return;
    }
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step("<<a<<","<<b<<"), N(b)="<<b.nm<<endl;
#endif
  Quad q = a/b;  // rounded, so N(a/b -q) is minimal
#ifdef DEBUG_PSEA
  cout<<" - translation = "<<q<<endl;
#endif
  Quad u;
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
      cout << endl;
#endif
      if (compute_t) t=0;
      return;
    }

  // Now look or a suitable alpha, trying all in turn (skipping alpha=0)

  Quad r,s, a1,b1, bs;
  int local_t = 1;
  for (vector<mat22>::iterator Mi=M_alphas.begin()+1; Mi!=M_alphas.end(); Mi++, local_t++)
    {
#ifdef DEBUG_PSEA
      cout<<" - testting type "<<t<<", M="<<(*Mi)<<endl;
#endif
      mat22 M = *Mi;
      r=-M.d, s=M.c; // alpha = r/s
      // Find the shift taking a/b closest to alpha
      bs = b*s;
      q = (a*s-b*r)/bs; // closest integer to (a/b)-(r/s)
      // We need to use temporary copies of a,b in case this alpha fails
      a1 = a-q*b, b1 = b;
      Mi->apply_left(a1,b1);
      if (b1.nm < b.nm) // success!
        {
          a = a1;
          b = b1;
          if (compute_c1d1)
            {
              d1 += q*c1;
              Mi->apply_right_inverse(c1,d1);
            }
          if (compute_c2d2)
            {
              d2 += q*c2;
              Mi->apply_right_inverse(c2,d2);
            }
#ifdef DEBUG_PSEA
          cout<<" - success, returning (a,b) = ("<<a<<","<<b<<")"<<endl;
#endif
          if (compute_t) t=local_t;
          return;
        }
#ifdef DEBUG_PSEA
      else
        {
          cout<<" - failure, new b would have had norm "<<quadnorm(b1)<<endl;
        }
#endif
    }
  // We should never arrive here as it means that all alphas have failed
  if (compute_t) t = -1;
  cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
  exit(1);
}

Quad quadgcd_psea(const Quad& aa, const Quad& bb)   // Using (pseudo-)EA
{
  Quad a(aa), b(bb);
  while (b.nm) pseudo_euclidean_step(a, b);
   while (!pos(a)) a*=fundunit;
 return a;
}

Quad quadbezout_psea(const Quad& aa, const Quad& bb, Quad& xx, Quad& yy)   // Using (pseudo-)EA
{
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
