// FILE EUCLID.CC: implementation of pseudo-euclidean algorithm and associated data

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

  // if (d==67)
  //   {
  //     Quad den(3,-1);
  //     alphas.push_back(RatQuad(6+w,den));
  //     alphas.push_back(RatQuad(-6-w,den));
  //     alphas.push_back(RatQuad(2+w,den));
  //     alphas.push_back(RatQuad(-2-w,den));
  //     den = quadconj(den);
  //     alpha_denoms.push_back(den);
  //     alphas.push_back(RatQuad(7-w,den));
  //     alphas.push_back(RatQuad(w-7,den));
  //     alphas.push_back(RatQuad(3-w,den));
  //     alphas.push_back(RatQuad(w-3,den));

  //     alphas.push_back(RatQuad(w,4));
  //     alphas.push_back(RatQuad(-w,4));
  //     alphas.push_back(RatQuad(w-1,4));
  //     alphas.push_back(RatQuad(1-w,4));
  //     alphas.push_back(RatQuad(1+w,4));
  //     alphas.push_back(RatQuad(-1-w,4));
  //     alphas.push_back(RatQuad(w-2,4));
  //     alphas.push_back(RatQuad(2-w,4));
  //     n_alphas += 16;
  //     return;
  //   }
  cerr << "define_alphas() not yet implemented for field "<<d<<endl;
}


// pseudo-Euclidean step: applies a translation and M_alpha inversion
// to a/b (or column vector [a;b]) reducing b, also multiplying row
// vector [c.d] my M_alpha on the right.  In the Euclidean case, the
// shift is -q where q=a/b (rounded) and the inversion is via
// S=[0,-1;1,0].

// a,b,c,d are changed in place, and on return, t holds the "type"
// (index of alpha which worked)

//#define DEBUG_PSEA

void pseudo_euclidean_step(Quad& a, Quad& b, Quad& c, Quad& d, int& t)
{
  t = 0;
  long normb = quadnorm(b);
  if (normb==0)
    return;
#ifdef DEBUG_PSEA
  cout<<"Entering pseudo_euclidean_step("<<a<<","<<b<<"), N(b)="<<normb<<endl;
#endif
  Quad q = a/b;  // rounded, so N(a/b -q) is minimal
#ifdef DEBUG_PSEA
  cout<<" - translation = "<<q<<endl;
#endif
  Quad u;
  a -= q*b;
  d += q*c;
#ifdef DEBUG_PSEA
  cout<<" - reduced a = "<<a<<endl;
#endif
  if (quadnorm(a) < normb) // always true in Euclidean case; invert using S
    {
      u = a; a=-b; b=u;
      u = d; d=-c; c=u;
#ifdef DEBUG_PSEA
      cout<<" - now N(a)<N(b), returning (a,b) = ("<<a<<","<<b<<")"<<endl;
#endif
      return;
    }

  // Now look or a suitable alpha, trying all in turn (skipping alpha=0)

  Quad r,s, a1,b1, bs;
  t = 1;
  for (vector<mat22>::iterator Mi=M_alphas.begin()+1; Mi!=M_alphas.end(); Mi++, t++)
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
      if (quadnorm(b1) < normb) // success!
        {
          a = a1;
          b = b1;
          d += q*c;
          Mi->apply_right_inverse(c,d);
#ifdef DEBUG_PSEA
      cout<<" - success, returning (a,b) = ("<<a<<","<<b<<")"<<endl;
#endif
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
  t = -1;
  cerr<<"Pseudo-Euclidean step fails for ("<<a<<", "<<b<<")"<<endl;
  exit(1);
}

