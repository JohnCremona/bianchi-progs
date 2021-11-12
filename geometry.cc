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
// alpha_inv is a permutation of range(n_alphas) such that
// alpha_inv[i]=j where M_alpha[i](oo) = alpha[j].
//
// alpha_flip is a permutation of range(n_alphas) such that
// alpha_flip[i]=j where -alpha[i] = alpha[j] mod 1.
//
// sigma_flip is a permutation of range(n_sigmas) such that
// sigma_flip[i]=j where -sigma[i] = sigma[j] mod 1.


int n_alphas, n_sigmas;
vector<RatQuad> alphas;
vector<RatQuad> sigmas;
vector<mat22> M_alphas;
vector<int> alpha_inv;
vector<int> alpha_flip;
vector<int> sigma_flip;
vector<int> edge_pairs_minus;
vector<int> edge_pairs_plus;
vector<int> edge_fours;
vector<int> cyclic_triangles;
vector<vector<int> > triangles;
vector<pair<vector<int>, Quad>> aas_triangles;
vector<pair<vector<int>, vector<Quad>> > squares;
vector<pair<vector<int>, vector<Quad>> > hexagons;

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  RatQuad alpha(-d,c);
  mat22 M(a,b,c,d);  // maps alpha = -d/c to oo
  //cout<<"M_alpha = "<<M<<" with determinant "<<M.det()<<endl;
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

// add r/s to the list of singular points, and optionally also -r/s (unless -r/s=r/s mod 1)

void add_sigma(const Quad& r, const Quad& s, int both_signs=1)
{
  RatQuad sigma(r,s);
  sigmas.push_back(sigma);
  if (both_signs)
    {
      sigmas.push_back(-sigma);
      sigma_flip.push_back(n_sigmas+1);   // transposition with next
      sigma_flip.push_back(n_sigmas);     // transposition with previous
      n_sigmas+=2;
    }
  else
    {
      sigma_flip.push_back(n_sigmas);     // identity
      n_sigmas+=1;
    }
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
      Quad u = d/2 -1-w;
      add_alpha(1+w,u,2,-1-w);  // alpha[1] = (1+w)/2
      alpha_inv.push_back(1);   // 1-1
      alpha_flip.push_back(1);  // 1-1
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
        }
      break;
    }
  case 2: case 5:
    {
      if (d!=5)
        {
          add_alpha_pair(3, w, +1);
          add_sigma(w,3);
        }
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
      if (d>6)
        {
          add_alpha_foursome(3, w+1, w-1);
          add_sigma(w,3);
          break;
        }
    }
  } // d%12

  switch (d) {
  case 19:
  case 43:
  default:
    return;
  case 5:
    {
      add_alpha_pair(2*w, w-4, +1);      // N(s)=20
      return;
    }
  case 6:
    {
      add_alpha_pair(2*w, 5, +1);      // N(s)=24
      return;
    }
  case 23:
    {
      add_alpha_pair(w+1, 2-w, +1); // N(s)=8
      add_alpha_pair(2-w, 1+w, +1);
      add_alpha_pair(2+w, w-3, +1); // N(s)=12
      add_alpha_pair(3-w, -2-w, +1);
      return;
    }
  case 31:
    {
      add_alpha_pair(w, 3, +1);       // N(s)=8
      add_alpha_pair(1-w, 3, +1);
      add_alpha_pair(1+w, 3);         // N(s)=10
      add_alpha_pair(2-w, 3);
      add_alpha_pair(3+w, w-6, +1);   // N(s)=20
      add_alpha_pair(4-w, 5+w, +1);
      return;
    }
  case 67:
    {
      add_alpha_foursome(4, w, w-1);   // alphas  9,10,11,12
      add_alpha_foursome(4, w+1, 2-w); // alphas 13,14,15,16
      add_alpha_foursome(3-w, w+6, 2+w); // alphas 17,18,19,20
      add_alpha_foursome(2+w, 7-w, 3-w); // alphas 21,22,23,24
      assert (n_alphas==25);
      assert (M_alphas.size()==25);
      assert (alpha_inv.size()==25);
      return;
    }
  case 163:
    {
      add_alpha_foursome(4, w, w-1);   // alphas  9,10,11,12
      add_alpha_foursome(4, w+1, 2-w); // alphas 13,14,15,16
      assert (n_alphas==17);

      // For d=163 we have 82 more alphas!

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

    } // end of d=163 block
  }
}

/****************************************************************

 Code defining triangle relations

***************************************************************/

#define CHECK_TRIANGLES

int check_triangle(const vector<int>& T)
{
  // cout<<"Checking aaa-triangle "<<T<<endl;
  mat22 Mi=M_alphas[T[0]], Mj=M_alphas[T[1]], Mk=M_alphas[T[2]];
  RatQuad x = (Mi(Mj.preimage_oo()) - Mk.preimage_oo());
  return (x.is_integral());
}

int check_cyclic_triangle(int i)
{
  // cout<<"Checking cyclic triangle {"<<i<<"}"<<endl;
  Quad t=M_alphas[i].trace();
  return (t*t==1);
}

int check_aas_triangle(const vector<int>& T, const Quad& u)
{
  // cout<<"Checking aas-triangle ("<<T<<","<<u<<")"<<endl;
  int i=T[0], j=T[1], k=T[2];
  RatQuad x = M_alphas[i](sigmas[j]+u) - sigmas[k];
  return (x.is_integral());
}

// aaa triangle relations:

// {i,j,k} where M_alphas[i](alphas[j]) = x + alphas[k] with x integral.

// The triangle has vertices [alpha_i, oo, alpha_j] and edges
// (I)_i = {alpha_i, oo},
// (M1)_j' = M_j' * {alpha_j',oo} = {oo, alpha_j},
// (M2)_k = M_i' * T^x * {alpha_k, oo} = M_i' * {x+alpha_k, oo} = {alpha_j, alpha_i}

// aas triangle relations:

// ({i,j,k},u) where M_alphas[i](sigmas[j]+u) = x + sigmas[k] with x integral.

// The triangle has vertices [alpha_i, oo, sigma_j+u] and edges
// (I)_i = {alpha_i, oo},
// - (T^u)_{-j} = - T^u * {sigma_j,oo} = {oo, sigma_j+u},
// (M_i'*T^x)_{-k} = M_i' * T^x * {sigma_k, oo} = M_i' * {x+sigma_k, oo} = {sigma_j+u, alpha_i}

void make_triangles()
{
  // Lists of general triangles and cyclic triangles (aaa, i.e. all vertices principal cusps)
  // NB Not including the {0,1,oo} triangle which is handles as a special case.
  switch(Quad::d) {
  case 43:
    triangles = {{3,7,4}};
    cyclic_triangles = {};
    break;
  case 67:
    triangles = {{3,7,4}, {0,19,24}, {1,22,17}, {3,13,18}, {3,22,15}, {9,13,16}};
    cyclic_triangles = {9};
    break;
  case 163:
    triangles = {{3,7,4}, {3, 21, 49}, {3, 53, 23}, {3, 71, 87}, {3, 83, 94}, {3, 85, 79}, {3, 98, 84},
                 {7, 98, 93}, {9, 13, 16}, {9, 69, 78}, {13, 25, 77}, {13, 33, 57}, {13, 61, 30}, {13, 69, 26},
                 {17, 33, 29}, {21, 28, 32}, {21, 35, 27}, {21, 59, 79}, {21, 71, 63}, {21, 75, 94}, {21, 98, 67},
                 {23, 51, 3}, {25, 33, 23}, {25, 45, 67}, {27, 76, 48}, {27, 79, 13}, {29, 48, 93}, {31, 35, 19},
                 {33, 45, 98}, {37, 40, 44}, {41, 45, 48}, {41, 73, 66}, {47, 96, 33}, {49, 84, 67}, {49, 87, 63},
                 {53, 83, 75}, {53, 85, 59}, {59, 89, 62}, {61, 69, 23}, {73, 92, 21}, {77, 87, 5} };
    cyclic_triangles = {9, 17};
    break;
  case 23:
    triangles = {{0,3,6}, {0,5,8}, {0,9,12}, {1,6,12}, {7,9,2}};
    cyclic_triangles = {};
    break;
  case 31:
    triangles = {{0,3,10}, {0,5,12}, {0,7,14}, {0,15,7}, {1,19,10}, {3,6,8}, {3,10,13}, {3,15,12},
                 {7,11,10}, {7,19,17}};
    cyclic_triangles = {};
    break;
  case 5:
  case 6:
  default:
    triangles = {};
    cyclic_triangles = {};
  }

#ifdef CHECK_TRIANGLES
  for (vector<vector<int> >::const_iterator Ti = triangles.begin(); Ti!=triangles.end(); ++Ti)
    {
      assert(check_triangle(*Ti));
    }
  for (vector<int>::const_iterator Ti = cyclic_triangles.begin(); Ti!=cyclic_triangles.end(); ++Ti)
    assert(check_cyclic_triangle(*Ti));
#endif

  if (Quad::class_number==1)
    {
      aas_triangles = {};
      return;
    }

  // Lists of aas triangles (i.e. two principal and one non-principal vertex)
  Quad w = Quad::w;
  switch(Quad::d) {
  case 5:
    aas_triangles = {{{1,1,1},0}, {{2,1,1},0}};
    break;
  case 6:
    aas_triangles = {{{1,1,1},0}, {{3,1,1},0}};
    break;
  case 23:
    aas_triangles = {{{1,1,1},0}, {{6,1,1},0}, {{12,1,1},0}, {{1,2,2},w}, {{8,2,2},0}, {{10,2,2},0}};
    break;
  case 31:
    aas_triangles = {{{11,1,1},0}, {{2,1,1},-w}, {{2,2,2},-1}, {{10,2,2},w-1} };
    break;
  default:
    aas_triangles = {};
  }
#ifdef CHECK_TRIANGLES
  for (vector<pair<vector<int>, Quad>>::const_iterator Ti = aas_triangles.begin(); Ti!=aas_triangles.end(); ++Ti)
    {
      assert(check_aas_triangle(Ti->first, Ti->second));
    }
#endif
}


/****************************************************************

 Code defining square relations

***************************************************************/

#define CHECK_SQUARES

int check_square(const vector<int>& S, const vector<Quad>& xyz)
{
  // Check:  the square has vertices {alpha_i, oo, alpha[j']+z, beta}
  // where beta = z + M_j(x+alpha[k']) = M_i'(y+alpha_l),
  // so that M_i(T^z(M_j(x+alpha[k']))) = y+alpha_l.

  // Edges:

  // {alpha_i, oo} = (I)_i
  // {oo, alpha_j'+z} = (T^z*M_j)_j
  // {alpha_j'+z, beta} = (T^z*M_j*T^x*M_k)_k
  // {beta, alpha_i} = (M_i'*T^y)_l

  // int i=S[0], j=S[1], k=S[2], l=S[3];
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  mat22 Mi=M_alphas[S[0]], Mj=M_alphas[S[1]], Mk=M_alphas[S[2]], Ml=M_alphas[S[3]];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));  // = x+alpha_k'
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0)); // = y+alpha_l
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  return ((M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1)) == alpha2);
}

void make_squares()
{
  Quad w = Quad::w;
  // Lists of general squares
  switch(Quad::d) {
  case 19:
    squares = {{{0,1,0,1}, {0,0,0}}};
    break;
  case 43:
    squares = {{{0,5,1,6}, {1-w, 0,0}},
               {{1,7,1,7}, {1, -1,0}}};
    break;
  case 67:
    squares = {{{3, 2, 4,12}, {0,0,0}},
               {{9, 0, 9, 0}, {0,0,0}},
               {{13, 0,14, 8}, {0,1,0}}};
    break;
  case 163:
    squares = {{{42,36,8,30}, {0,0,0}},
               {{17,43,17,43}, {0,0,0}},
               {{38,0,37,19}, {0,0,0}},
               {{41,35,7,29}, {0,0,0}},
               {{11,31,89,33}, {0,0,0}},
               {{1,90,1,90}, {-w,w,w}},
               {{88,9,87,8}, {0,0,0}},
               {{1,89,1,89}, {1,-1,0}},
               {{62,8,59,2}, {-1,0,0}},
               {{11,17,11,17}, {0,0,0}},
               {{2,66,0,75}, {0,0,0}},
               {{50,9,55,2}, {0,0,0}},
               {{83,11,81,0}, {0,0,0}}};
    break;
  case 5:
    squares = {{{0, 1, 0, 1}, {0,0,0}},
               {{2, 1, 3, 1}, {-1,-1-w,0}},
               {{0, 3, 0, 2}, {-1,-1,0}}};
    break;
  case 6:
    squares = {{{0, 3, 0, 2}, {0,0,0}},
               {{3, 1, 2, 1}, {1,-w,0}}};
    break;
  default:
    squares = {};
  }

#ifdef CHECK_SQUARES
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Si = squares.begin(); Si!=squares.end(); ++Si)
    {
      // cout<<"Checking square "<<Si->first<<", "<<Si->second<<endl;
      assert(check_square(Si->first, Si->second));
    }
#endif
}

/****************************************************************

 Code defining hexagon relations

***************************************************************/

#define CHECK_HEXAGONS

int check_hexagon(const vector<int>& ijklmn, const vector<Quad>& ux1y1x2y2)
{
  // Check:  the hexagon has vertices {beta_1, alpha_i, oo, u+alpha[j], beta_2, gamma}
  // where beta1 = M_i'(x1+alpha[k]), beta2 = M_j'(x2+alpha[l]),
  // gamma = M_i'*T^x1*M_k'*T^y1(alpha[m]) = T^u*M_j'*T^x2*M_l'*T^y2(alpha[n]).

  // Edges:

  // +{alpha_i, oo} = (I)_i
  // +{beta1, alpha_i} = (M_i'*T^x1)_k
  // +{gamma, beta1} = (M_i'*T^x1*M_k'*T^y1)_m
  // -{u+alpha_j, oo} = - (T^u)_j
  // -{beta2, alpha_j} = - (T^u*M_j'*T^x2)_l
  // -{gamma, beta2} = - (T^u*M_j'*T^x2*M_l'*T^y2)_n

  Quad u = ux1y1x2y2[0], x1 = ux1y1x2y2[1], y1 = ux1y1x2y2[2], x2 = ux1y1x2y2[3], y2 = ux1y1x2y2[4];
  int i=ijklmn[0], j=ijklmn[1], k=ijklmn[2], l=ijklmn[3], m=ijklmn[4], n=ijklmn[5];
  RatQuad gamma1 = (M_alphas[alpha_inv[i]]*mat22::Tmat(x1)*M_alphas[alpha_inv[k]])(y1+alphas[m]);
  RatQuad gamma2 = (mat22::Tmat(u)*M_alphas[alpha_inv[j]]*mat22::Tmat(x2)*M_alphas[alpha_inv[l]])(y2+alphas[n]);
  return gamma1==gamma2;
}

void make_hexagons()
{
  Quad w = Quad::w;
  // Lists of general hexagons
  switch(Quad::d) {
  case 6:
    hexagons = {{{1, 0, 0, 1, 1, 0}, {w+1, w, -w-1, -w, w+1}}};
    break;
  default:
    hexagons = {};
  }

#ifdef CHECK_HEXAGONS
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Hi = hexagons.begin(); Hi!=hexagons.end(); ++Hi)
    {
      // cout<<"Checking hexagon "<<Hi->first<<", "<<Hi->second<<endl;
      assert(check_hexagon(Hi->first, Hi->second));
    }
#endif
}

void make_faces()
{
  make_triangles();
  make_squares();
  make_hexagons();
}

