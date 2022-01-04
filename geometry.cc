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
vector<TRIANGLE> aaa_triangles;
vector<TRIANGLE> aas_triangles;
vector<POLYGON> squares;
vector<POLYGON> hexagons;

// Global functions to be used once during setting the field:

// read from data file 'geodata.dat' into global variables alphas, sigmas,
// aaa_triangles, aas_triangles, cyclic_triangles, squares, hexagons
void read_data(int verbose=0);

void check_triangles(int verbose=0);
void check_squares(int verbose=0);
void check_hexagons(int verbose=0);
void alphas_sigmas_universal();
void alphas_sigmas_denom_2();
void alphas_sigmas_denom_3();

void Quad::setup_geometry()
{
  n_alphas = n_sigmas = 0;

  alphas_sigmas_universal();

  if (Quad::is_Euclidean) return;

  alphas_sigmas_denom_2();
  alphas_sigmas_denom_3();

  int debug = 0;
  read_data(debug); // read remaining alphas and sigmas and all faces from geodat.dat
  if (debug)
    {
      cout << "alphas:\n";
      for(int i=0; i<n_alphas; i++) cout<<i<<": "<<alphas[i]<<endl;
      cout << "sigmas:\n";
      for(int i=0; i<n_sigmas; i++) cout<<i<<": "<<sigmas[i]<<endl;
    }
#ifdef CHECK_TRIANGLES
  check_triangles(debug);
#endif
#ifdef CHECK_SQUARES
  check_squares(debug);
#endif
#ifdef CHECK_HEXAGONS
  check_hexagons(debug);
#endif
}

// Given a,b,c,d with ad-bc=2 add alpha=-d/c and M_alpha=[a,b;c,d] to the global lists

void add_alpha(const Quad& a, const Quad& b, const Quad& c, const Quad& d)
{
  RatQuad alpha(-d,c);
  mat22 M(a,b,c,d);  // maps alpha = -d/c to oo
  //cout<<"M_alpha = "<<M<<" with determinant "<<M.det()<<endl;
  assert (M.is_unimodular());
  alphas.push_back(alpha);
  M_alphas.push_back(M);
  n_alphas++;
}

// If r1*r2 = -1 mod s with r1, r2 distinct we have r1/s, -r1/s, r2/s,
// -r2/s with matrices [r2,t;s,-r1] and similar.  When r1==r2 we have
// aa "minus pair" with r1^2=-1 (mod s), and when r1==-r2 we have a
// "plus pair" with r1^2=+1 (mod s).

void add_alpha_orbit(const Quad& s, const Quad& r1, const Quad& r2)
{
  assert(div(s,r1*r2+Quad::one));
  Quad t = -(r1*r2+Quad::one)/s;
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

void add_sigma_orbit(const Quad& r, const Quad& s)
{
  RatQuad sigma(r,s);
  sigmas.push_back(sigma);
  if (s.is_zero() || s==Quad(2)) // don't also include -sigma
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

void alphas_sigmas_universal()
{
  // sigma_0 = oo:

  add_sigma_orbit(Quad::one,Quad::zero); // fill in the 0'th entry in sigmas,
                  // so the others will be indexed from 1

  // alpha_0 = 0:

  add_alpha(Quad::zero,-Quad::one,Quad::one,Quad::zero);  // alpha[0] = 0
  alpha_inv.push_back(0); // 0-0
  alpha_flip.push_back(0); // 0-0
  assert (n_alphas==1);
}

// alphas and sigmas with denominator 2:

void alphas_sigmas_denom_2()
{
  long d = Quad::d;
  QUINT D(d);
  Quad w = Quad::w, two(2);

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
      Quad u = (D-1)/2;
      add_alpha(w,u,two,-w);  // alpha[1] = w/2
      alpha_inv.push_back(1); // 1-1
      alpha_flip.push_back(1); // 1-1
      add_sigma_orbit(w+Quad::one,two);
      break;
    }
  case 2: // (2) = (2,w)^2
  case 6:
    {
      Quad u = D/2 -1-w;
      add_alpha(Quad::one+w,u,two,-Quad::one-w);  // alpha[1] = (1+w)/2
      alpha_inv.push_back(1);   // 1-1
      alpha_flip.push_back(1);  // 1-1
      add_sigma_orbit(w,two);
      break;
    }
  case 3:
    {
      Quad u = (D-3)/8;  // = 2, 5, 8, 20 for d=19,43,67,163 = 3 (mod 8) so 2 is inert
      add_alpha(w-Quad::one,u,two,-w);  // alpha[1] = w/2
      add_alpha(w,u,two,Quad::one-w);   // alpha[2] = (w-1)/2
      alpha_inv.push_back(2); // 1-2
      alpha_inv.push_back(1); // 2-1
      alpha_flip.push_back(1); // 1-1
      alpha_flip.push_back(2); // 2-2
      break;
    }
  case 7: // (2) = (2,w)*(2,1-w)
    {
      add_sigma_orbit(w,two);
      add_sigma_orbit(Quad::one-w,two);
      break;
    }
  } // d%8
}

// alphas and sigmas with denominator 3:

void alphas_sigmas_denom_3()
{
  int d = Quad::d;
  Quad w = Quad::w, three(3);

  switch (d%12) {
  case 1: case 10:
    {
      add_alpha_orbit(three, w, w);
      add_alpha_orbit(three, Quad::one+w, Quad::one-w);
      // no sigmas (3 inert)
      break;
    }
  case 7:
    {
      if (d>19) // e.g. 43, 67, 163
        {
          add_alpha_orbit(three, w, Quad::one-w);
          add_alpha_orbit(three, Quad::one+w, Quad::one+w);
        }
      // no sigmas (3 inert)
      break;
    }
  case 2: case 5:
    {
      if (d!=5)
        {
          add_alpha_orbit(three, w, -w);
          add_sigma_orbit(Quad::one+w,three);
          add_sigma_orbit(Quad::one-w,three);
        }
      break;
    }
  case 11:
    {
      add_alpha_orbit(three, Quad::one+w, -Quad::one-w);
      if (d==35)
        {
          add_sigma_orbit(w,three);
        }
      if (d>35)
        {
          add_sigma_orbit(w,three);
          add_sigma_orbit(w-Quad::one,three);
        }
      break;
    }
  case 3:
    {
      if (d>15)
        {
          add_alpha_orbit(three, w, w-Quad::one);
          add_sigma_orbit(Quad::one+w,three);
        }
      break;
    }
  case 6:  case 9:
    {
      if (d>6)
        {
          add_alpha_orbit(three, w+Quad::one, w-Quad::one);
          add_sigma_orbit(w,three);
          break;
        }
    }
  } // d%12
}


/****************************************************************

 Code defining triangle relations

***************************************************************/

int check_aaa_triangle(const TRIANGLE& T, int verbose)
{
  const vector<int>& tri = T.first;
  const Quad& u = T.second;
  if (verbose)
    cout<<"Checking aaa-triangle ("<<tri<<","<<u<<")"<<endl;
  mat22 Mi=M_alphas[tri[0]], Mj=M_alphas[tri[1]], Mk=M_alphas[tri[2]];
  RatQuad x = (Mi(Mj.preimage_oo()+u) - Mk.preimage_oo());
  if (verbose)
    cout<<"x = "<<x<<endl;
  return (x.is_integral());
}

int check_cyclic_triangle(int i, int verbose)
{
  if (verbose)
     cout<<"Checking cyclic triangle {"<<i<<"}"<<endl;
  Quad t=M_alphas[i].trace();
  return (t*t==Quad::one);
}

int check_aas_triangle(const TRIANGLE& T, int verbose)
{
  const vector<int>& tri = T.first;
  const Quad& u = T.second;
  if (verbose)
    cout<<"Checking aas-triangle ("<<tri<<","<<u<<")"<<endl;
  int i=tri[0], j=tri[1], k=tri[2];
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

void check_triangles(int verbose)
{
  for (vector<int>::const_iterator Ti = cyclic_triangles.begin(); Ti!=cyclic_triangles.end(); ++Ti)
    assert(check_cyclic_triangle(*Ti, verbose));
  for (vector<TRIANGLE>::const_iterator Ti = aaa_triangles.begin(); Ti!=aaa_triangles.end(); ++Ti)
    {
      if(verbose)
        cout<<"Checking aaa-triangle, ijk="<<Ti->first<<", u="<<Ti->second<<endl;
      assert(check_aaa_triangle(*Ti, verbose));
    }
  for (vector<TRIANGLE>::const_iterator Ti = aas_triangles.begin(); Ti!=aas_triangles.end(); ++Ti)
    {
      if(verbose)
        cout<<"Checking aas-triangle, ijk="<<Ti->first<<", u="<<Ti->second<<endl;
      assert(check_aas_triangle(*Ti, verbose));
    }
}

/****************************************************************

 Code defining square relations

***************************************************************/

int check_square(const POLYGON& squ)
{
  // Check:  the square {{i,j,k,l},{x,y,z}} has vertices {alpha_i, oo, alpha[j']+z, beta}
  // where beta = z + M_j(x+alpha[k']) = M_i'(y+alpha_l),
  // so that M_i(T^z(M_j(x+alpha[k']))) = y+alpha_l.

  // Edges:

  // {alpha_i, oo} = (I)_i
  // {oo, alpha_j'+z} = (T^z*M_j)_j
  // {alpha_j'+z, beta} = (T^z*M_j*T^x*M_k)_k
  // {beta, alpha_i} = (M_i'*T^y)_l

  const vector<int>& S = squ.first;  // int i=S[0], j=S[1], k=S[2], l=S[3];
  const vector<Quad>& xyz = squ.second;
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  mat22 Mi=M_alphas[S[0]], Mj=M_alphas[S[1]], Mk=M_alphas[S[2]], Ml=M_alphas[S[3]];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));  // = x+alpha_k'
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0)); // = y+alpha_l
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  return ((M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1)) == alpha2);
}

void check_squares(int verbose)
{
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Si = squares.begin(); Si!=squares.end(); ++Si)
    {
      if (verbose)
        cout<<"Checking square "<<Si->first<<", "<<Si->second<<endl;
      assert(check_square(*Si));
    }
}

/****************************************************************

 Code defining hexagon relations

***************************************************************/

int check_hexagon(const POLYGON& hex)
{
  // Check: the hexagon {{i,j,k,l,m,n},{u,x1,y1,x2,y2}} has vertices
  // {beta_1, alpha_i, oo, u+alpha[j], beta_2, gamma} where

  // beta1 = M_i'(x1+alpha[k]),
  // beta2 = M_j'(x2+alpha[l]),
  // gamma = M_i'*T^x1*M_k'*T^y1(alpha[m]) = T^u*M_j'*T^x2*M_l'*T^y2(alpha[n]).

  // Edges:

  // +{alpha_i, oo} = (I)_i
  // +{beta1, alpha_i} = (M_i'*T^x1)_k
  // +{gamma, beta1} = (M_i'*T^x1*M_k'*T^y1)_m
  // -{u+alpha_j, oo} = - (T^u)_j
  // -{beta2, alpha_j} = - (T^u*M_j'*T^x2)_l
  // -{gamma, beta2} = - (T^u*M_j'*T^x2*M_l'*T^y2)_n

  const vector<int>& ijklmn = hex.first;
  const vector<Quad>& ux1y1x2y2 = hex.second;
  Quad u = ux1y1x2y2[0], x1 = ux1y1x2y2[1], y1 = ux1y1x2y2[2], x2 = ux1y1x2y2[3], y2 = ux1y1x2y2[4];
  int i=ijklmn[0], j=ijklmn[1], k=ijklmn[2], l=ijklmn[3], m=ijklmn[4], n=ijklmn[5];
  RatQuad gamma1 = (M_alphas[alpha_inv[i]]*mat22::Tmat(x1)*M_alphas[alpha_inv[k]])(y1+alphas[m]);
  RatQuad gamma2 = (mat22::Tmat(u)*M_alphas[alpha_inv[j]]*mat22::Tmat(x2)*M_alphas[alpha_inv[l]])(y2+alphas[n]);
  return gamma1==gamma2;
}

void check_hexagons(int verbose)
{
  for (vector<pair<vector<int>, vector<Quad>> >::const_iterator Hi = hexagons.begin(); Hi!=hexagons.end(); ++Hi)
    {
      if (verbose)
        cout<<"Checking hexagon "<<Hi->first<<", "<<Hi->second<<endl;
      assert(check_hexagon(*Hi));
    }
}

/****************************************************************

 Code for reading face relations from file 'geodata.dat'

***************************************************************/

// reads from data file into global variables aaa_triangles,
// aas_triangles, cyclic_triangles, squares, hexagons

// Each line starts with d, with d==0 meaning "ignore this line" (for
// comments) and d==-1 means "stop reading", with the d values in
// increasing order so that we can stop reading when we see a d value
// larger than Quad::d.

// The next non-space character is A,S,T,U,C,Q,H:

// A 'alpha orbit' followed by 6 integers: s, r1, r2
// S 'sigma orbit' followed by 4 integers: r, s
// T 'aaa-triangle' followed by 5 integers: i,j,k; u
// U 'aas-triangle' followed by 5 integers: i,j,k; u
// C 'cyclic triangle' followed by 1 integer: i
// Q 'square' followed by 10 integers: i,j,k,l; x,y,z
// H 'hexagon' followed by 16 integers: i,j,k,l,m,n; u,x1,y1,x2,y2

void read_data(int verbose)
{
  ifstream geodata;
  string line;
  int file_d, i, j, k, l, m, n;
  char G;
  Quad s, r, r1, r2, u, x, y, z, x1, y1, x2, y2;

  geodata.open("geodata.dat");
  getline(geodata, line);
  while (!geodata.eof())
    {
      istringstream input_line(line);
      input_line >> file_d;
      if (file_d==-1 || file_d > Quad::d)
        break;
      if (file_d != Quad::d)
        {
          if (verbose>1)
            cout<<"Skipping line: "<<line<<endl;
          getline(geodata, line);
          continue;
        }
      if (verbose)
        cout<<"Processing line "<<line<<endl;
      input_line >> G;
      switch(G) {
      case 'A': // alpha orbit
        {
          if (verbose)
            cout << " reading alpha orbit data"<<endl;
          input_line >> s >> r1 >> r2;
          if (verbose)
            cout << " - (s,r1,r2) = ("<< s<<","<<r1<<","<<r2<<")" <<endl;
          add_alpha_orbit(s, r1, r2);
          break;
        }
      case 'S': // sigma orbit
        {
          if (verbose)
            cout << " reading sigma orbit data"<<endl;
          input_line >> r >> s;
          if (verbose)
            cout << " - (r,s) = ("<< r <<","<< s <<")" <<endl;
          add_sigma_orbit(r, s);
          break;
        }
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
      case 'T': // aaa-triangle
      case 'U': // aas-triangle
        {
          if (verbose)
            cout << " reading AA"<<(G=='T'? 'A': 'S')<<"-triangle data"<<endl;
          input_line >> i >> j >> k >> u;
          if (verbose)
            cout << " - [i,j,k] = ["<<i<<","<<j<<","<<k<<"], u = "<<u<<endl;
          if (G=='T')
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

