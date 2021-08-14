// FILE GEOMETRY.CC: implementation of functions and associated data for hyperbolic tessation

#include <iostream>

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
// alpha_flip is a permutation of range(N_alphas) such that
// alpha_flip[i]=j where -alpha[i] = alpha[j] mod 1.


int n_alphas, n_sigmas;
vector<RatQuad> alphas;
vector<RatQuad> sigmas;
vector<mat22> M_alphas;
vector<int> alpha_inv;
vector<int> alpha_flip;
vector<int> edge_pairs_minus;
vector<int> edge_pairs_plus;
vector<int> edge_fours;
vector<int> cyclic_triangles;
vector<vector<int> > triangles;
vector<pair<vector<int>, Quad>> aas_triangles;
vector<pair<vector<int>, vector<Quad>> > squares;

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  RatQuad alpha(-d,c);
  mat22 M(a,b,c,d);  // maps alpha = -d/c to oo
  //  cout<<"M_alpha = "<<M<<" with determinant "<<M.det()<<endl;
  assert (M.det()==1);
  alphas.push_back(alpha);
  M_alphas.push_back(M);
  n_alphas++;
}

// If r1*r2 = -1 mod s with r1, r2 distinct we have r1/s, -r1/s, r2/s, -r2/s
// with matrices [r2,t;s,-r1] and similar

void add_alpha_foursome(const Quad& s, const Quad& r1, const Quad& r2)
{
  edge_fours.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+2);
  alpha_inv.push_back(n_alphas+3);
  alpha_inv.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+1);
  alpha_flip.push_back(n_alphas+1);
  alpha_flip.push_back(n_alphas);
  alpha_flip.push_back(n_alphas+3);
  alpha_flip.push_back(n_alphas+2);
  Quad t = -(r1*r2+1)/s;
  add_alpha( r2, t, s, -r1); // alpha =  r1/s
  add_alpha(-r2, t, s,  r1); // alpha = -r1/s
  add_alpha( r1, t, s, -r2); // alpha =  r2/s
  add_alpha(-r1, t, s,  r2); // alpha = -r2/s
}

// If r*r = -1 mod s we have r/s, -r/s
// with matrices [r,t;s,-r] and [-r,t;s,r]
// where t = (-r^2-1)/s.
// If r*r = +1 mod s we have r/s, -r/s
// with matrices [-r,t;s,-r] and [r,t;s,r]
// where t = (r^2-1)/s.

void add_alpha_pair(const Quad& s, const Quad& r, int sign=-1)
{
  if (sign==-1)
    {
      edge_pairs_minus.push_back(n_alphas);
      alpha_inv.push_back(n_alphas);   // identity
      alpha_inv.push_back(n_alphas+1); // identity
    }
  else
    {
      edge_pairs_plus.push_back(n_alphas);
      alpha_inv.push_back(n_alphas+1); // transposition with next
      alpha_inv.push_back(n_alphas);   // transposition with previous
    }
  alpha_flip.push_back(n_alphas+1);   // transposition with next
  alpha_flip.push_back(n_alphas);     // transposition with previous
  Quad t = (r*r*sign-1)/s;
  add_alpha(-sign*r, t, s, -r); // alpha =  r/s
  add_alpha( sign*r, t, s,  r); // alpha = -r/s
}

// add r/s to the list of singular points, and optionally also -r/s

void add_sigma(const Quad& r, const Quad& s, int both_signs=1)
{
  RatQuad sigma(r,s);
  sigmas.push_back(sigma);
  n_sigmas++;
  if (both_signs)
    {
      sigmas.push_back(-sigma);
      n_sigmas++;
    }
}

#define CHECK_TRIANGLES

int check_triangle(int i, int j, int k)
{
  mat22 Mi=M_alphas[i];
  RatQuad aj = M_alphas[j].inverse()(RatQuad(1,0)),
    ak = M_alphas[k].inverse()(RatQuad(1,0));
  RatQuad x = (Mi(aj) - ak);
  return (x.is_integral());
}

void add_triangle(int i, int j, int k)
{
  triangles.push_back({i,j,k});
#ifdef CHECK_TRIANGLES
  assert (check_triangle(i,j,k));
#endif
}

int check_cyclic_triangle(int i)
{
  Quad t=M_alphas[i].trace();
  return (t*t==1);
}

void add_cyclic_triangle(int i)
{
  cyclic_triangles.push_back(i);
#ifdef CHECK_TRIANGLES
  assert (check_cyclic_triangle(i));
#endif
}

int check_aas_triangle(int i, int j, int k, const Quad& u)
{
  RatQuad x = M_alphas[i](sigmas[j]+u) - sigmas[k];
  return (x.is_integral());
}

void add_aas_triangle(int i, int j, int k, const Quad& u=0)
{
  aas_triangles.push_back({{i,j,k},u});
#ifdef CHECK_TRIANGLES
  assert (check_aas_triangle(i, j, k, u));
#endif
}

#define CHECK_SQUARES

int check_square(int i, int j, int k, int l, const Quad& x, const Quad& y, const Quad& z)
{
  // Check:  the square has vertices {alpha_i, oo, alpha[j'], beta}
  // where beta = z + M_j(x+alpha[k']) = M_i'(y+alpha_l),
  // so that M_i(T^z(M_j(x+alpha[k']))) = y+alpha_l.

  // Edges:

  // {alpha_i, oo} = (I)_i
  // {oo, alpha_j'+z} = (T^z*M_j)_j
  // {alpha_j'+z, beta} = (T^z*M_j*T^k*M_k)_k
  // {beta, alpha_i} = (M_i'*T^y)_l

  mat22 Mi=M_alphas[i], Mj=M_alphas[j], Mk=M_alphas[k], Ml=M_alphas[l];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));  // = x+alpha_k'
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0)); // = y+alpha_l
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  return ((M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1)) == alpha2);
}

void add_square(int i, int j, int k, int l, const Quad& x=Quad::zero, const Quad& y=Quad::zero, const Quad& z=Quad::zero)
{
  vector<int> squ = {i,j,k,l};
  vector<Quad> xyz = {x,y,z};
  squares.push_back({squ,xyz});
#ifdef CHECK_SQUARES
  assert (check_square(i,j,k,l, x,y,z));
#endif
}

// Global function to be used once during setting the field:

void Quad::setup_geometry()
{
  int d = Quad::d;
  n_alphas = n_sigmas = 0;

  add_sigma(1,0, 0); // fill in the 0'th entry in sigmas,
                  // so the others will be indexed from 1

  // alphas (only 0) with denominator 1:

  add_alpha(0,-1,1,0);  // alpha[0] = 0
  alpha_inv.push_back(0); // 0-0
  alpha_flip.push_back(0); // 0-0
  assert (n_alphas==1);

  if (Quad::is_Euclidean) return;

  Quad w = Quad::w;

  // alphas and sigmas with denominator 2:

  // these are always w/2 and (1+w)/2 but we use (w-1)/2 instead of (w+1)/2 when d%4=3

  // alpha = w/2, sigma = (w+1)/2 when d%4=1, 2 ramifies, (2)=(2,1+w)^2
  // alpha = (w+1)/2, sigma = w/2 when d%4=2, 2 ramifies, (2)=(2,w)^2
  // alpha = w/2, (w-1)/2 when 2 is inert, d%8=3
  // sigma = w/2, (w-1)/2 when 2 splits, d%8=7

  // These are alphas number (up to) 1 and 2 and sigmas (up to) 1 and
  // 2, and determine edge orbit numbers (up to) 1 and 2, and (up to)
  // n_alphas and n_alphas+1 respectively.

  switch (d%8) {
  case 1: // (2) = (2,w+1)^2
  case 5:
    {
      Quad u = (d-1)/2;
      add_alpha(w,u,2,-w);  // alpha[1] = w/2
      alpha_inv.push_back(1); // 1-1
      alpha_flip.push_back(1); // 1-1
      add_sigma(w+1,2, 0);
      break;
    }
  case 2: // (2) = (2,w)^2
  case 6:
    {
      Quad u = d/2;
      add_alpha(1+w,u,2,-1-w);  // alpha[1] = (1+w)/2
      alpha_flip.push_back(1); // 1-1
      add_sigma(w,2, 0);
      break;
    }
  case 3:
    {
      Quad u = (d-3)/8;  // = 2, 5, 8, 20 for d=19,43,67,163 = 3 (mod 8) so 2 is inert
      add_alpha(w-1,u,2,-w);  // alpha[1] = w/2
      add_alpha(w,u,2,1-w);   // alpha[2] = (w-1)/2
      alpha_inv.push_back(2); // 1-2
      alpha_inv.push_back(1); // 2-1
      alpha_flip.push_back(1); // 1-1
      alpha_flip.push_back(2); // 2-2
      break;
    }
  case 7: // (2) = (2,w)*(2,1-w)
    {
      add_alpha_pair(4, 1+2*w, +1);
      add_sigma(w,2, 0);
      add_sigma(1-w,2, 0);
      break;
    }
  } // d%8

  // alphas and sigmas with denominator 3:

  switch (d%12) {
  case 1: case 10:
    {
      add_alpha_pair(3, w, -1);
      add_alpha_foursome(3, 1+w, 1-w);
      break;
    }
  case 7:
    {
      if (d>19) // e.g. 43, 67, 163
        {
          add_alpha_foursome(3, w, 1-w);
          add_alpha_pair(3, 1+w, -1);
          add_triangle(3,7,4); // <w/3, oo, (w+1)/3>
        }
      break;
    }
  case 2: case 5:
    {
      add_alpha_pair(3, w, +1);
      add_sigma(w,3);
      break;
    }
  case 11:
    {
      add_alpha_pair(3, 1+w, +1);
      if (d>23)
        {
          add_sigma(w,3);
          add_sigma(1-w,3);
        }
      break;
    }
  case 3:
    {
      add_alpha_foursome(3, w, w-1);
      add_sigma(1+w,3);
      break;
    }
  case 6:
  case 9:
    {
      add_alpha_foursome(3, w+1, w-1);
      add_sigma(w,3);
      break;
    }
  } // d%12

  if (d==5)
    {
      add_alpha_pair(2*w, 4+w, +1);      // N(s)=20
      return;
    }

  if (d==19) // no more alphas
    {
      add_square(0,1,0,1);  // symmetric
      return;
    }

  if (d==23)
    {
      add_alpha_pair(w+1, 2-w, +1); // N(s)=8
      add_alpha_pair(2-w, 1+w, +1);
      add_alpha_pair(2+w, w-3, +1); // N(s)=12
      add_alpha_pair(3-w, -2-w, +1);

      // AAA-triangles (all vertices principal "Alpha"s)
      add_triangle(0,3,6);
      add_triangle(0,5,8);
      add_triangle(0,9,12);
      add_triangle(1,6,12);
      add_triangle(7,9,2);
      // AAS-triangles (two vertices principal "Alpha"s, one non-principal "Sigma")
      add_aas_triangle(1,1,1);
      add_aas_triangle(6,1,1);
      add_aas_triangle(12,1,1);
      add_aas_triangle(1,2,2, w);
      add_aas_triangle(8,2,2);
      add_aas_triangle(10,2,2);

      // cout<<n_alphas<<" alphas: "<<alphas<<endl;
      // cout<<"alpha_inv: "<<alpha_inv<<endl;
      // cout<<"alpha_flip: "<<alpha_flip<<endl;
      // cout<<"aaa-triangles: ";
      // for (vector<vector<int>>::const_iterator Ti=triangles.begin(); Ti!=triangles.end(); Ti++)
      //   cout<<(*Ti)<<" ";
      // cout<<endl;
      // cout<<n_sigmas<<" sigmas: "<<sigmas<<endl;
      // cout<<"aas-triangles: ";
      // for (vector<pair<vector<int>, Quad>>::const_iterator Ti=aas_triangles.begin(); Ti!=aas_triangles.end(); Ti++)
      //   cout<<"["<<Ti->first<<"; "<<Ti->second<<"] ";
      // cout<<endl;
      return;
    }

  if (d==43)
    {
      add_square(0,5,1,6, 1-w,  0);
      add_square(1,7,1,7,   1, -1); // symmetric
      return;
    }

  if (d==31)
    {
      add_alpha_pair(w, 3, +1);       // N(s)=8
      add_alpha_pair(1-w, 3, +1);
      add_alpha_pair(1+w, 3);         // N(s)=10
      add_alpha_pair(2-w, 3);
      add_alpha_pair(3+w, w-6, +1);   // N(s)=20
      add_alpha_pair(4-w, 5+w, +1);
      return;
    }

  // alphas with denominator 4 for both d=67 and d=163:

  if ((d==67)||(d==163))
    {
      add_alpha_foursome(4, w, w-1);   // alphas  9,10,11,12
      add_alpha_foursome(4, w+1, 2-w); // alphas 13,14,15,16
      assert (n_alphas==17);
    }

  if (d==67)   // alphas with denominator norm 23 for d=67 only:
    {
      add_alpha_foursome(3-w, w+6, 2+w); // alphas 17,18,19,20
      add_alpha_foursome(2+w, 7-w, 3-w); // alphas 21,22,23,24
      assert (n_alphas==25);
      assert (M_alphas.size()==25);
      assert (alpha_inv.size()==25);

      add_cyclic_triangle(9);
      add_triangle(0,19,24);
      add_triangle(1,22,17);
      add_triangle(3,13,18);
      add_triangle(3,22,15);
      add_triangle(9,13,16);

      add_square( 3, 2, 4,12, 0);
      add_square( 9, 0, 9, 0, 0); // symmetric
      add_square(13, 0,14, 8, 0, 1);

      return;
    }

  if (d==163)
    {  // For d=163 we have 82 more alphas!

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

      // Add 42 extra triangles

      // cyclic triangles [9, 10, 11, 12, 17, 18, 19, 20]
      add_cyclic_triangle(9);
      add_cyclic_triangle(17);

      // add_triangle(3,7,4); // <w/3, oo, (w+1)/3> already added above
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

      add_square(42,36,8,30);
      add_square(17,43,17,43);     // symmetric
      add_square(38,0,37,19);
      add_square(41,35,7,29);
      add_square(11,31,89,33);
      add_square(1,90,1,90, -w,w,w);
      add_square(88,9,87,8);
      add_square(1,89,1,89, 1,-1); // symmetric
      add_square(62,8,59,2, -1,0);
      add_square(11,17,11,17);     // symmetric
      add_square(2,66,0,75);
      add_square(50,9,55,2);
      add_square(83,11,81,0);

    } // end of d=163 block
}
