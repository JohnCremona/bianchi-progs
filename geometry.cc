// FILE GEOMETRY.CC: implementation of functions and associated data for hyperbolic tessation

#include <iostream>

#include "geometry.h"

#define CHECK_TRIANGLES
#define CHECK_SQUARES
#define CHECK_HEXAGONS


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

// Definitions of global objects declared in geometry.h:

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
vector<pair<vector<int>, Quad>> aaa_triangles;
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

// If r1*r2 = -1 mod s with r1, r2 distinct we have r1/s, -r1/s, r2/s,
// -r2/s with matrices [r2,t;s,-r1] and similar.  When r1==r2 we have
// aa "minus pair" with r1^2=-1 (mod s), and when r1==-r2 we have a
// "plus pair" with r1^2=+1 (mod s).

void add_alpha_foursome(const Quad& s, const Quad& r1, const Quad& r2)
{
  Quad t = -(r1*r2+1)/s;
  if (r1==r2) // "-" pair, r1*r2=-1 (mod s)
    {
      edge_pairs_minus.push_back(n_alphas);
      alpha_inv.push_back(n_alphas);   // identity
      alpha_inv.push_back(n_alphas+1); // identity
      alpha_flip.push_back(n_alphas+1);   // transposition with next
      alpha_flip.push_back(n_alphas);     // transposition with previous
      add_alpha( r1, t, s, -r1); // alpha =  r1/s
      add_alpha(-r1, t, s,  r1); // alpha = -r1/s
      return;
    }
  if (r1==-r2) // "+" pair, r1*r2=+1 (mod s)
    {
      edge_pairs_plus.push_back(n_alphas);
      alpha_inv.push_back(n_alphas+1); // transposition with next
      alpha_inv.push_back(n_alphas);   // transposition with previous
      alpha_flip.push_back(n_alphas+1);   // transposition with next
      alpha_flip.push_back(n_alphas);     // transposition with previous
      add_alpha(-r1, t, s, -r1); // alpha =  r1/s
      add_alpha( r1, t, s,  r1); // alpha = -r1/s
      return;
    }
  // Now we have four distinct alphas

  edge_fours.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+2);
  alpha_inv.push_back(n_alphas+3);
  alpha_inv.push_back(n_alphas);
  alpha_inv.push_back(n_alphas+1);
  alpha_flip.push_back(n_alphas+1);
  alpha_flip.push_back(n_alphas);
  alpha_flip.push_back(n_alphas+3);
  alpha_flip.push_back(n_alphas+2);
  add_alpha( r2, t, s, -r1); // alpha =  r1/s
  add_alpha(-r2, t, s,  r1); // alpha = -r1/s
  add_alpha( r1, t, s, -r2); // alpha =  r2/s
  add_alpha(-r1, t, s,  r2); // alpha = -r2/s
}

// add r/s to the list of singular points, and optionally also -r/s (unless -r/s=r/s mod 1)

void add_sigma(const Quad& r, const Quad& s)
{
  RatQuad sigma(r,s);
  sigmas.push_back(sigma);
  if (s==0 || s==2) // don't also include -sigma
    {
      sigma_flip.push_back(n_sigmas);     // identity
      n_sigmas+=1;
    }
  else
    {
      sigmas.push_back(-sigma);
      sigma_flip.push_back(n_sigmas+1);   // transposition with next
      sigma_flip.push_back(n_sigmas);     // transposition with previous
      n_sigmas+=2;
    }
}

// Global functions to be used once during setting the field:

// read from data file 'geodata.dat' into global variables alphas, sigmas,
// aaa_triangles, aas_triangles, cyclic_triangles, squares, hexagons
void read_data(int verbose=0);

void check_triangles();
void check_squares();
void check_hexagons();
void alphas_sigmas_universal();
void alphas_sigmas_denom_2();
void alphas_sigmas_denom_3();
void alphas_sigmas_other();

void Quad::setup_geometry()
{
  n_alphas = n_sigmas = 0;

  alphas_sigmas_universal();

  if (Quad::is_Euclidean) return;

  alphas_sigmas_denom_2();
  alphas_sigmas_denom_3();
  alphas_sigmas_other();

  read_data(0);
#ifdef CHECK_TRIANGLES
  check_triangles();
#endif
#ifdef CHECK_SQUARES
  check_squares();
#endif
#ifdef CHECK_HEXAGONS
  check_hexagons();
#endif
}

void alphas_sigmas_universal()
{
  // sigma_0 = oo:

  add_sigma(1,0); // fill in the 0'th entry in sigmas,
                  // so the others will be indexed from 1

  // alpha_0 = 0:

  add_alpha(0,-1,1,0);  // alpha[0] = 0
  alpha_inv.push_back(0); // 0-0
  alpha_flip.push_back(0); // 0-0
  assert (n_alphas==1);
}

// alphas and sigmas with denominator 2:

void alphas_sigmas_denom_2()
{
  int d = Quad::d;
  Quad w = Quad::w;

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
      add_sigma(w+1,2);
      break;
    }
  case 2: // (2) = (2,w)^2
  case 6:
    {
      Quad u = d/2 -1-w;
      add_alpha(1+w,u,2,-1-w);  // alpha[1] = (1+w)/2
      alpha_inv.push_back(1);   // 1-1
      alpha_flip.push_back(1);  // 1-1
      add_sigma(w,2);
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
      add_alpha_foursome(4, 1+2*w, -(1+2*w));
      add_sigma(w,2);
      add_sigma(1-w,2);
      break;
    }
  } // d%8
}

// alphas and sigmas with denominator 3:

void alphas_sigmas_denom_3()
{
  int d = Quad::d;
  Quad w = Quad::w;

  switch (d%12) {
  case 1: case 10:
    {
      add_alpha_foursome(3, w, w);
      add_alpha_foursome(3, 1+w, 1-w);
      break;
    }
  case 7:
    {
      if (d>19) // e.g. 43, 67, 163
        {
          add_alpha_foursome(3, w, 1-w);
          add_alpha_foursome(3, 1+w, 1+w);
        }
      break;
    }
  case 2: case 5:
    {
      if (d!=5)
        {
          add_alpha_foursome(3, w, -w);
          add_sigma(w,3);
        }
      break;
    }
  case 11:
    {
      add_alpha_foursome(3, 1+w, -1-w);
      if (d>23)
        {
          add_sigma(w,3);
          add_sigma(w-1,3);
        }
      break;
    }
  case 3:
    {
      add_alpha_foursome(3, w, w-1);
      add_sigma(1+w,3);
      break;
    }
  case 6:  case 9:
    {
      if (d>6)
        {
          add_alpha_foursome(3, w+1, w-1);
          add_sigma(w,3);
          break;
        }
    }
  } // d%12
}

void alphas_sigmas_other()
{
  int d = Quad::d;
  Quad w = Quad::w;

  switch (d) {
  case 19:
  case 43:
  default:
    return;
  case 5:
    {
      add_alpha_foursome(2*w, w-4, 4-w);      // N(s)=20
      return;
    }
  case 6:
    {
      add_alpha_foursome(2*w, 5, -5);      // N(s)=24
      return;
    }
  case 10:
    {
      add_alpha_foursome(2*w, 9, -9);      // N(s)=40
      return;
    }
  // case 13:
  //   {
  //     add_alpha_foursome(5, 2*w + 2, 2*w - 2);
  //     return;
  //   }
  // case 14:
  //   {
  //     add_alpha_foursome(w - 2, w + 5, w + 3);
  //     add_alpha_foursome(w + 2, w - 3, w - 5);
  //     return;
  //   }
  // case 15:
  //   {
  //     return;
  //   }
  // case 17:
  //   {
  //     add_alpha_foursome(w - 1, 7, 5);
  //     add_alpha_foursome(w + 1, 7, 5);
  //     add_alpha_foursome(7, 3*w - 3, 3*w + 3);
  //     return;
  //   }
  // case 21:
  //   {
  //     add_alpha_foursome(w - 1, 9, -5);
  //     add_alpha_foursome(w + 1, 9, -5);
  //     add_alpha_foursome(9, 4*w + 4, 4*w - 4);
  //     return;
  //   }
  // case 22:
  //   {
  //     add_alpha_foursome(4, w + 1, -w + 1);
  //     add_alpha_foursome(w, 9, -5);
  //     add_alpha_foursome(w - 1, 9, 5);
  //     add_alpha_foursome(w + 1, 9, 5);
  //     add_alpha_foursome(5, w + 2, w - 2);
  //     add_alpha_foursome(5, 2*w + 1, -2*w + 1);
  //     add_alpha_foursome(w + 2, w - 5, w - 9);
  //     add_alpha_foursome(w - 2, w + 9, w + 5);
  //     return;
  //   }
  case 23:
    {
      add_alpha_foursome(w+1, 2-w, w-2); // N(s)=8
      add_alpha_foursome(2-w, 1+w, -1-w);
      add_alpha_foursome(2+w, w-3, 3-w); // N(s)=12
      add_alpha_foursome(3-w, -2-w, w+2);
      return;
    }
  case 31:
    {
      add_alpha_foursome(w, 3, -3);       // N(s)=8
      add_alpha_foursome(1-w, 3, -3);
      add_alpha_foursome(1+w, 3, 3);         // N(s)=10
      add_alpha_foursome(2-w, 3, 3);
      add_alpha_foursome(3+w, w-6, 6-w);   // N(s)=20
      add_alpha_foursome(4-w, 5+w, -5-w);
      return;
    }
  case 47:
    {
      add_alpha_foursome(w - 1, 5, -5);
      add_alpha_foursome(w, 5, -5);
      add_alpha_foursome(w - 4, w + 3, -3-w);
      add_alpha_foursome(w - 4, 2*w + 3, -3-2*w);
      add_alpha_foursome(w + 3, w - 4, 4-w);
      add_alpha_foursome(w + 3, 2*w - 5, 5-2*w);
      add_alpha_foursome(6, 2*w - 1, 1-2*w);
      add_alpha_foursome(w - 3, w + 4, w + 2);
      add_alpha_foursome(w + 2, w - 3, w - 5);
      add_alpha_foursome(5, 2*w + 1, 2*w - 3);
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
      add_alpha_foursome(7, 2+3*w, 2+3*w);
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

int check_aaa_triangle(const vector<int>& T, const Quad& u)
{
  // cout<<"Checking aaa-triangle ("<<T<<","<<u<<")"<<endl;
  mat22 Mi=M_alphas[T[0]], Mj=M_alphas[T[1]], Mk=M_alphas[T[2]];
  RatQuad x = (Mi(Mj.preimage_oo()+u) - Mk.preimage_oo());
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

// {{i,j,k},u} where M_alphas[i](alphas[j]+u) = x + alphas[k] with x integral.

// The triangle has vertices [alpha_i, oo, alpha_j+u] and edges
// (I)_i = {alpha_i, oo},
// (M1)_j' = T^u * M_j' * {alpha_j',oo} = {oo, alpha_j +u},
// (M2)_k = M_i' * T^x * {alpha_k, oo} = M_i' * {x+alpha_k, oo} = {alpha_j +u, alpha_i}

// aas triangle relations:

// ({i,j,k},u) where M_alphas[i](sigmas[j]+u) = x + sigmas[k] with x integral.

// The triangle has vertices [alpha_i, oo, sigma_j+u] and edges
// (I)_i = {alpha_i, oo},
// - (T^u)_{-j} = - T^u * {sigma_j,oo} = {oo, sigma_j+u},
// (M_i'*T^x)_{-k} = M_i' * T^x * {sigma_k, oo} = M_i' * {x+sigma_k, oo} = {sigma_j+u, alpha_i}

void check_triangles()
{
  for (vector<pair<vector<int>, Quad>>::const_iterator Ti = aaa_triangles.begin(); Ti!=aaa_triangles.end(); ++Ti)
    {
      assert(check_aaa_triangle(Ti->first, Ti->second));
    }
  for (vector<int>::const_iterator Ti = cyclic_triangles.begin(); Ti!=cyclic_triangles.end(); ++Ti)
    assert(check_cyclic_triangle(*Ti));
  for (vector<pair<vector<int>, Quad>>::const_iterator Ti = aas_triangles.begin(); Ti!=aas_triangles.end(); ++Ti)
    {
      assert(check_aas_triangle(Ti->first, Ti->second));
    }
}

/****************************************************************

 Code defining square relations

***************************************************************/

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

void check_squares()
{
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Si = squares.begin(); Si!=squares.end(); ++Si)
    {
      // cout<<"Checking square "<<Si->first<<", "<<Si->second<<endl;
      assert(check_square(Si->first, Si->second));
    }
}

/****************************************************************

 Code defining hexagon relations

***************************************************************/

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

void check_hexagons()
{
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Hi = hexagons.begin(); Hi!=hexagons.end(); ++Hi)
    {
      // cout<<"Checking hexagon "<<Hi->first<<", "<<Hi->second<<endl;
      assert(check_hexagon(Hi->first, Hi->second));
    }
}

/****************************************************************

 Code for reading face relations from file 'geodata.dat'

***************************************************************/

// reads from data file into global variables aaa_triangles,
// aas_triangles, cyclic_triangles, squares, hexagons

void read_data(int verbose)
{
  ifstream geodata;
  string line;
  int file_d, i, j, k, l, m, n;
  char G;
  Quad u, x, y, z, x1, y1, x2, y2;

  geodata.open("geodata.dat");
  getline(geodata, line);
  while (!geodata.eof())
    {
      istringstream input_line(line);
      input_line >> file_d;
      if (file_d==-1)
        break;
      if (file_d != Quad::d)
        {
          if (verbose)
            cout<<"Skipping line: "<<line<<endl;
          getline(geodata, line);
          continue;
        }
      if (verbose)
        cout<<"Processing line "<<line<<endl;
      input_line >> G;
      switch(G) {
      case 'C': // cyclic triangle
        {
          if (verbose)
            cout << " reading cyclic triangle data"<<endl;
          input_line >> j;
          if (verbose)
            cout << " - j = "<< j <<endl;
          cyclic_triangles.push_back(j);
          break;
        }
      case 'A': // aaa-triangle
      case 'S': // aas-triangle
        {
          if (verbose)
            cout << " reading AA"<<G<<"-triangle data"<<endl;
          input_line >> i >> j >> k >> u;
          if (verbose)
            cout << " - [i,j,k] = ["<<i<<","<<j<<","<<k<<"], u = "<<u<<endl;
          if (G=='A')
            aaa_triangles.push_back({{i,j,k}, u});
          else
            aas_triangles.push_back({{i,j,k}, u});
          break;
        }
      case 'Q': // square
        {
          if (verbose)
            cout << " reading square data"<<endl;
          input_line >> i >> j >> k >> l >> x >> y >> z;
          if (verbose)
            cout << " - [i,j,k,l] = ["<<i<<","<<j<<","<<k<<","<<l<<"], "
                 << "[x,y,z] = ["<<x<<","<<y<<","<<z<<"]"<<endl;
          squares.push_back({{i,j,k,l}, {x,y,z}});
          break;
        }
      case 'H': // hexagon
        {
          if (verbose)
            cout << " reading hexagon data"<<endl;
          input_line >> i >> j >> k >> l >> m >> n >> u >> x1 >> y1 >> x2 >> y2;
          if (verbose)
            cout << " - [i,j,k,l,m,n] = ["<<i<<","<<j<<","<<k<<","<<l<<","<<m<<","<<n<<"], "
                 << "[u,x1,y1,x2,y2] = ["<<u<<","<<x1<<","<<y1<<","<<x2<<","<<y2<<"]"<<endl;
          hexagons.push_back({{i,j,k,l,m,n}, {u,x1,y1,x2,y2}});
          break;
        }
      default:
        break;
      }
      getline(geodata, line);
    }
  geodata.close();
}

