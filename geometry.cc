// FILE GEOMETRY.CC: implementation of functions and associated data for hyperbolic tessation

#include <iostream>

#include <eclib/arith.h>
#include "quads.h"
#include "ratquads.h"
#include "geometry.h"

// Definitions of commonly used matrices

mat22 mat22::identity(1,0,0,1);
mat22 mat22::J(-1,0,0,1);
mat22 mat22::S(0,-1,1,0);
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
vector<int> edge_pairs;
vector<int> edge_fours;
vector<int> cyclic_triangles;
vector<vector<int> > triangles;
vector<pair<vector<int>, vector<Quad>> > squares;

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

void add_edge_foursome(const Quad& s, const Quad& r1, const Quad& r2)
{
  edge_fours.push_back(n_alphas);
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
  edge_pairs.push_back(n_alphas);
  Quad t = -(r*r+1)/s;
  alpha_inv.push_back(n_alphas);
  add_alpha( r, t, s, -r); // alpha =  r/s
  alpha_inv.push_back(n_alphas);
  add_alpha(-r, t, s,  r); // alpha = -r/s
}

void add_triangle(int i, int j, int k)
{
  triangles.push_back({i,j,k});
  // Check
  mat22 Mi=M_alphas[i], Mj=M_alphas[j], Mk=M_alphas[k];
  RatQuad beta(-Mj.entry(1,1),Mj.entry(1,0));
  RatQuad gamma(-Mk.entry(1,1),Mk.entry(1,0));
  RatQuad x = (Mi.entry(0,0)*beta+Mi.entry(0,1))/(Mi.entry(1,0)*beta+Mi.entry(1,1)) - gamma;
  assert (x.is_integral());
}

void add_cyclic_triangle(int i)
{
  cyclic_triangles.push_back(i);
  // Check
  Quad t=M_alphas[i].trace();
  assert (t*t==1);
}

void add_square(int i, int j, int k, int l, const Quad& x, const Quad& y)
{
  vector<int> squ = {i,j,k,l};
  vector<Quad> xy = {x,y};
  squares.push_back({squ,xy});
  // Check
  mat22 Mi=M_alphas[i], Mj=M_alphas[j], Mk=M_alphas[k], Ml=M_alphas[l];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0));
  mat22 M = Mi*Mj;
  assert ((M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1)) == alpha2);
}

// Global function to be used once during setting the field:

void Quad::setup_geometry()
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

  if (d==19)
    {
      add_square(0,1,0,1, 0, 0);
      return;
    }

  // alphas (w/3, -w/3, (1-w)/3, (w-1)/3, (w+1)/3, -(w+1)/3) with denominator 3:

  add_edge_foursome(3, w, 1-w);
  add_alpha_pair(3, 1+w);
  assert (n_alphas==9);
  assert (M_alphas.size()==9);
  assert (alpha_inv.size()==9);

  add_triangle(3,7,4); // <w/3, oo, (w+1)/3>

  if (d==43)
    {
      add_square(0,5,1,6, 1-w,  0);
      add_square(1,7,1,7,   1, -1);
      return;
    }

  // alphas with denominator 4 for both d=67 and d=163:

  add_edge_foursome(4, w, w-1);   // alphas  9,10,11,12
  add_edge_foursome(4, w+1, 2-w); // alphas 13,14,15,16
  assert (n_alphas==17);

  if (d==67)   // alphas with denominator norm 23 for d=67 only:
    {
      add_edge_foursome(3-w, w+6, 2+w); // alphas 17,18,19,20
      add_edge_foursome(2+w, 7-w, 3-w); // alphas 21,22,23,24
      assert (n_alphas==25);
      assert (M_alphas.size()==25);
      assert (alpha_inv.size()==25);

      add_cyclic_triangle(9);
      add_triangle(0,19,24);
      add_triangle(1,22,17);
      add_triangle(3,13,18);
      add_triangle(3,22,15);
      add_triangle(9,13,16);

      add_square( 3, 2, 4,12, 0, 0);
      add_square( 9, 0, 9, 0, 0, 0);
      add_square(13, 0,14, 8, 0, 1);

      return;
    }

  // Now d=163 and we have 82 more alphas!

  // 20 alphas with s=5 (norm 25)

  add_edge_foursome(5, w, w-1);
  add_edge_foursome(5, 2*w, 2-2*w);
  add_edge_foursome(5, w+2, 1-2*w);
  add_edge_foursome(5, w-2, 2+2*w);
  add_edge_foursome(5, w+1, 1+2*w);
  assert (n_alphas==37);

  // 12 alphas with s=6 (norm 36)

  add_edge_foursome(6, w, 1-w);
  add_edge_foursome(6, 1+w, w-2);
  add_edge_foursome(6, 2+w, 3-w);
  assert (n_alphas==49);

  // 4 alphas with s=w (norm 41), and 4 conjugates of these

  add_edge_foursome(w, 12, 17);
  add_edge_foursome(1-w, 12, 17);
  assert (n_alphas==57);

  // 4 alphas with s=1+w (norm 43), and 4 conjugates of these

  add_edge_foursome(1+w, 12, w-17);
  add_edge_foursome(2-w, 12, -w-16);
  assert (n_alphas==65);

  // 8 alphas with s=2+w (norm 47), and 8 conjugates of these

  add_edge_foursome(2+w, 7, 18-w);
  add_edge_foursome(2+w, w-11, w-16);

  add_edge_foursome(3-w, 7, 17+w);
  add_edge_foursome(3-w, -w-10, -w-15);
  assert (n_alphas==81);

  // 10 alphas with s=7 (norm 49)

  add_edge_foursome(7, w+3, 2*w-1);
  add_edge_foursome(7, 2*w+1, 3-2*w);
  add_alpha_pair(7, 2+3*w);
  assert (n_alphas==91);

  // 4 alphas with s=3+w (norm 53), and 4 conjugates of these

  add_edge_foursome(3+w, 5-w, w-17);
  add_edge_foursome(4-w, w+4, -w-16);
  assert (n_alphas==99);
  assert (M_alphas.size()==99);
  assert (alpha_inv.size()==99);

  // Add 42 extra triangles

  add_cyclic_triangle(9);
  add_cyclic_triangle(17);
  add_triangle(3, 21, 49);
  add_triangle(3, 53, 23);
  add_triangle(3, 71, 87);
  add_triangle(3, 83, 94);
  add_triangle(3, 85, 79);
  add_triangle(3, 98, 84);
  add_triangle(7, 98, 93);
  add_triangle(9, 13, 16);
  add_triangle(9, 69, 78);
  add_triangle(13, 25, 77);
  add_triangle(13, 33, 57);
  add_triangle(13, 61, 30);
  add_triangle(13, 69, 26);
  add_triangle(17, 33, 29);
  add_triangle(21, 28, 32);
  add_triangle(21, 35, 27);
  add_triangle(21, 59, 79);
  add_triangle(21, 71, 63);
  add_triangle(21, 75, 94);
  add_triangle(21, 98, 67);
  add_triangle(23, 51, 3);
  add_triangle(25, 33, 23);
  add_triangle(25, 45, 67);
  add_triangle(27, 76, 48);
  add_triangle(27, 79, 13);
  add_triangle(29, 48, 93);
  add_triangle(31, 35, 19);
  add_triangle(33, 45, 98);
  add_triangle(37, 40, 44);
  add_triangle(41, 45, 48);
  add_triangle(41, 73, 66);
  add_triangle(47, 96, 33);
  add_triangle(49, 84, 67);
  add_triangle(49, 87, 63);
  add_triangle(53, 83, 75);
  add_triangle(53, 85, 59);
  add_triangle(59, 89, 62);
  add_triangle(61, 69, 23);
  add_triangle(73, 92, 21);
  add_triangle(77, 87, 5);

}
