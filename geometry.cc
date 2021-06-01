// FILE GEOMETRY.CC: implementation of functions and associated data for hyperbolic tessation

#include <iostream>

#include <eclib/arith.h>
#include "quads.h"
#include "geometry.h"

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
// alpha_inv is a permutation of range(N_alphas) such that
// alpha_inv[i]=j where M_alpha[i](oo) = alpha[j].
//
// We do not store the alphas explicitly, they are -d/c where M_alpha=[a,b;c,d].

int n_alphas;
vector<mat22> M_alphas;
vector<int> alpha_inv;
vector<int> alpha_pairs;
vector<int> alpha_fours;

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  mat22 M_alpha(a,b,c,d);  // maps alpha = -d/c to oo
  assert (M_alpha.det()==1);
  M_alphas.push_back(M_alpha);
  n_alphas += 1;
}

// If r1*r2 = -1 mod s with r1, r2 distinct we have r1/s, -r1/s, r2/s, -r2/s
// with matrices [r2,t;s,-r1] and similar

void add_alpha_foursome(const Quad& s, const Quad& r1, const Quad& r2)
{
  alpha_fours.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+2);
  alpha_inv.push_back(n_alphas+3);
  alpha_inv.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+1);
  Quad t = -(r1*r2+1)/s;
  add_alpha( r2, t, s, -r1); // alpha =  r1/s
  add_alpha(-r2, t, s,  r1); // alpha = -r1/s
  add_alpha( r1, t, s, -r2); // alpha =  r2/s
  add_alpha(-r1, t, s,  r2); // alpha = -r2/s
}

// If r*r = -1 mod s we have r/s, -r/s
// with matrices [r,t;s,-r] and [-r,t;s,r]

void add_alpha_pair(const Quad& s, const Quad& r)
{
  alpha_pairs.push_back(n_alphas);
  Quad t = -(r*r+1)/s;
  alpha_inv.push_back(n_alphas);
  add_alpha( r, t, s, -r); // alpha =  r/s
  alpha_inv.push_back(n_alphas);
  add_alpha(-r, t, s,  r); // alpha = -r/s
}

// Global function to be used once after setting the field:

void define_alphas()
{
  int d = Quad::d;

  // alphas (only 0) with denominator 1:

  add_alpha(0,-1,1,0);  // alpha[0] = 0
  alpha_inv.push_back(0); // 0-0
  assert (n_alphas==1);

  if (d<19) return;

  Quad w = Quad::w;

  // alphas (w/2, (w-1)/2) with denominator 2:

  Quad u = (d-3)/8;  // = 2, 5, 8, 20
  add_alpha(w-1,u,2,-w);  // alpha[1] = w/2
  add_alpha(w,u,2,1-w);   // alpha[2] = (w-1)/2
  alpha_inv.push_back(2); // 1-2
  alpha_inv.push_back(1); // 2-1
  assert (n_alphas==3);
  assert (M_alphas.size()==3);
  assert (alpha_inv.size()==3);

  if (d<43) return;

  // alphas (w/3, -w/3, (1-w)/3, (w-1)/3, (w+1)/3, -(w+1)/3) with denominator 3:

  add_alpha_foursome(3, w, 1-w);
  add_alpha_pair(3, 1+w);
  assert (n_alphas==9);
  assert (M_alphas.size()==9);
  assert (alpha_inv.size()==9);

  if (d<67) return;

  // alphas with denominator 4 for both d=67 and d=163:

  add_alpha_foursome(4, w, w-1);
  add_alpha_foursome(4, w+1, 2-w);
  assert (n_alphas==17);

  if (d==67)   // alphas with denominator norm 23 for d=67 only:
    {
      add_alpha_foursome(3-w, w+6, 2+w);
      add_alpha_foursome(2+w, 7-w, 3-w);
      assert (n_alphas==25);
      assert (M_alphas.size()==25);
      assert (alpha_inv.size()==25);
      return;
    }

  // Now d=163 and we have 82 more alphas!

  // 20 alphas with s=5 (norm 25)

  add_alpha_foursome(5, w, w-1);
  add_alpha_foursome(5, 2*w, 2-2*w);
  add_alpha_foursome(5, w+2, 1-2*w);
  add_alpha_foursome(5, w-2, 2+2*w);
  add_alpha_foursome(5, w+1, 1+2*w);
  assert (n_alphas==37);

  // 12 alphas with s=6 (norm 36)

  add_alpha_foursome(6, w, 1-w);
  add_alpha_foursome(6, 1+w, w-2);
  add_alpha_foursome(6, 2+w, 3-w);
  assert (n_alphas==49);

  // 4 alphas with s=w (norm 41), and 4 conjugates of these

  add_alpha_foursome(w, 12, 17);
  add_alpha_foursome(1-w, 12, 17);
  assert (n_alphas==57);

  // 4 alphas with s=1+w (norm 43), and 4 conjugates of these

  add_alpha_foursome(1+w, 12, w-17);
  add_alpha_foursome(2-w, 12, -w-16);
  assert (n_alphas==65);

  // 8 alphas with s=2+w (norm 47), and 8 conjugates of these

  add_alpha_foursome(2+w, 7, 18-w);
  add_alpha_foursome(2+w, w-11, w-16);

  add_alpha_foursome(3-w, 7, 17+w);
  add_alpha_foursome(3-w, -w-10, -w-15);
  assert (n_alphas==81);

  // 10 alphas with s=7 (norm 49)

  add_alpha_foursome(7, w+3, 2*w-1);
  add_alpha_foursome(7, 2*w+1, 3-2*w);
  add_alpha_pair(7, 2+3*w);
  assert (n_alphas==91);

  // 4 alphas with s=3+w (norm 53), and 4 conjugates of these

  add_alpha_foursome(3+w, 5-w, w-17);
  add_alpha_foursome(4-w, w+4, -w-16);
  assert (n_alphas==99);
  assert (M_alphas.size()==99);
  assert (alpha_inv.size()==99);

}
