// FILE GEOMETRY.CC: implementation of functions and associated data for hyperbolic tessation

#include <iostream>

#include "geometry.h"
#include "swan.h"

SwanData Quad::SD;

/****************************************************************

 Code for reading face relations from file 'geodata.dat'

***************************************************************/
void Quad::setup_geometry(int debug)
{
  if (geometry_initialised)
    return;
  Quad::SD.read_geodata(0, "geodata", debug);
}

// reads from data file into global variables aaa_triangles,
// aas_triangles, squares, hexagons

// Each line starts with d, with d==0 meaning "ignore this line" (for
// comments) and d==-1 means "stop reading", with the d values in
// increasing order so that we can stop reading when we see a d value
// larger than Quad::d.

// The next non-space character is A,S,T,U,Q,H:

// A 'alpha orbit' followed by 6 integers: s, r1, r2
// S 'sigma orbit' followed by 4 integers: r, s
// T 'aaa-triangle' followed by 5 integers: i,j,k; u
// U 'aas-triangle' followed by 5 integers: i,j,k; u
// Q 'square' followed by 10 integers: i,j,k,l; x,y,z
// H 'hexagon' followed by 16 integers: i,j,k,l,m,n; u,x1,y1,x2,y2

void parse_geodata_line(const string& line, int& file_d, char& G, POLYGON& poly, int verbose)
{
  istringstream input_line(line);
  input_line >> file_d;
  if (file_d==-1 || file_d > Quad::d)
    {
      G = 'X'; // code for "quit reading"
      return;
    }
  if (file_d != Quad::d)
    {
      G = '%'; // code for "skip line"
      return;
    }
  if (verbose)
    cout<<"Processing line "<<line<<endl;
  input_line >> G;
  switch(G) {
  case 'A': // alpha orbit
    {
      Quad s, r1, r2;
      input_line >> s >> r1 >> r2;
      if (verbose)
        cout << " reading alpha orbit data"<<endl
             << " - (s,r1,r2) = ("<< s<<","<<r1<<","<<r2<<")" <<endl;
      poly = {{},{s,r1,r2}};
      return;
    }
  case 'S': // sigma orbit
    {
      Quad r, s;
      input_line >> r >> s;
      if (verbose)
        cout << " reading sigma orbit data"<<endl
             << " - (r,s) = ("<< r <<","<< s <<")" <<endl;
      poly = {{},{r,s}};
      return;
    }
  case 'T': // aaa-triangle
  case 'U': // aas-triangle
    {
      int i, j, k;   Quad u;
      input_line >> i >> j >> k >> u;
      if (verbose)
        cout << " reading AA"<<(G=='T'? 'A': 'S')<<"-triangle data"<<endl
             << " - [i,j,k] = ["<<i<<","<<j<<","<<k<<"], u = "<<u<<endl;
      poly = {{i,j,k}, {u}};
      return;
    }
  case 'Q': // square
    {
      int i, j, k, l; Quad x, y, z;
      input_line >> i >> j >> k >> l >> x >> y >> z;
      if (verbose)
        cout << " reading square data"<<endl
             << " - [i,j,k,l] = ["<<i<<","<<j<<","<<k<<","<<l<<"], "
             << "[x,y,z] = ["<<x<<","<<y<<","<<z<<"]"<<endl;
      poly = {{i,j,k,l}, {x,y,z}};
      return;
    }
  case 'H': // hexagon
    {
      int i, j, k, l, m, n; Quad u, x1, x2, y1, y2;
      input_line >> i >> j >> k >> l >> m >> n >> u >> x1 >> y1 >> x2 >> y2;
      if (verbose)
        cout << " reading hexagon data"<<endl
             << " - [i,j,k,l,m,n] = ["<<i<<","<<j<<","<<k<<","<<l<<","<<m<<","<<n<<"], "
             << "[u,x1,y1,x2,y2] = ["<<u<<","<<x1<<","<<y1<<","<<x2<<","<<y2<<"]"<<endl;
      poly = {{i,j,k,l,m,n}, {u,x1,y1,x2,y2}};
      return;
    }
  default:
    return;
  }
}

// same as above but only reads the polygons (T,U,Q,H):
// returns 4 lists, of aaa-triangles, aas-triangles, squares, hexagons
vector<vector<POLYGON>> read_polygons(string subdir, int verbose)
{
  ifstream geodata;
  string line;
  int file_d;
  char G;
  vector<POLYGON> Ts, Us, Qs, Hs;
  POLYGON poly;

  stringstream ss;
  if (subdir.empty())
    subdir = "geodata";
  ss << subdir << "/geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str());
  if (!geodata.is_open())
    {
      ss.clear();
      ss << "geodata_" << Quad::d << ".dat";
      geodata.open(ss.str().c_str());
      if (!geodata.is_open())
        {
          cout << "No geodata file!" <<endl;
          return {};
        }
    }
  if (verbose)
    cout << "reading from " << ss.str() <<endl;

  getline(geodata, line);
  while (!geodata.eof())
    {
      parse_geodata_line(line, file_d, G, poly, verbose);
      istringstream input_line(line);
      switch(G) {
      case 'X':
        {
          if (verbose>1)
            cout<<"Skipping rest of file"<<endl;
          break;
        }
      case '%':
      case 'A': // alpha orbit
      case 'S': // sigma orbit
      default:
        {
          if (verbose>1)
            cout<<"Skipping line: "<<line<<endl;
          break;
        }
      case 'T': // aaa-triangle
        {
          Ts.push_back(poly);
          break;
        }
      case 'U': // aas-triangle
        {
          Us.push_back(poly);
          break;
        }
      case 'Q': // square
        {
          Qs.push_back(poly);
          break;
        }
      case 'H': // hexagon
        {
          Hs.push_back(poly);
          break;
        }
      } // end of switch
      if (G=='X')
        break; // don't read any more lines
      else
        getline(geodata, line);
    }
  geodata.close();
  return {Ts, Us, Qs, Hs};
}

/***************************** obsolete code using globals ***************/
#if(0)

#define CHECK_TRIANGLES
#define CHECK_SQUARES
#define CHECK_HEXAGONS

// Definitions of alphas and associated matrices M_alpha such that
// det(M_alpha)=1 and M_alpha(alpha)=oo and M_alpha(oo)=alpha' in alphas.
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
map<vector<RAT>, int> alpha_ind;
map<vector<RAT>, int> sigma_ind;
std::set<Quad> alpha_denoms;
vector<mat22> M_alphas;
vector<int> alpha_inv;
vector<int> alpha_flip;
vector<int> sigma_flip;
vector<int> edge_pairs_minus;
vector<int> edge_pairs_plus;
vector<int> edge_fours;
vector<POLYGON> aaa_triangles;
vector<POLYGON> aas_triangles;
vector<POLYGON> squares;
vector<POLYGON> hexagons;

// Global functions to be used once during setting the field:

// read from data file 'geodata/geodata_{d}.dat' into global variables alphas, sigmas,
// aaa_triangles, aas_triangles, squares, hexagons
void read_data(int verbose=0);

int check_triangles(int verbose=0);
int check_squares(int verbose=0);
int check_hexagons(int verbose=0);
void alphas_sigmas_universal();
void alphas_sigmas_denom_2();
void alphas_sigmas_denom_3();

void Quad::setup_geometry(int debug)
{
  if (geometry_initialised)
    return;

  geometry_initialised = 1;
  n_alphas = n_sigmas = 0;

  // Clear these in case this is not the first field run

  alphas.clear();
  sigmas.clear();
  alpha_denoms.clear();
  M_alphas.clear();
  alpha_inv.clear();
  alpha_ind.clear();
  sigma_ind.clear();
  alpha_flip.clear();
  sigma_flip.clear();
  edge_pairs_plus.clear();
  edge_pairs_minus.clear();
  edge_fours.clear();
  aaa_triangles.clear();
  aas_triangles.clear();
  squares.clear();
  hexagons.clear();

  alphas_sigmas_universal();

  if (Quad::is_Euclidean) return;

  alphas_sigmas_denom_2();
  alphas_sigmas_denom_3();

  if (debug)
    {
      cout << "alphas and sigmas with denominators 1,2,3:" << endl;
      cout << "alphas:\n";
      for(int i=0; i<n_alphas; i++) cout<<i<<": "<<alphas[i]<<endl;
      cout << "sigmas:\n";
      for(int i=0; i<n_sigmas; i++) cout<<i<<": "<<sigmas[i]<<endl;
    }
  read_data(debug); // read remaining alphas and sigmas and all faces from geodat.dat
  if (debug)
    {
      cout << "alpha_denoms:\n";
      for (auto s : alpha_denoms)
        cout << s << " with norm "<<s.norm()<<endl;
      cout << "alphas:\n";
      for(int i=0; i<n_alphas; i++)
        cout<<i<<": "<<alphas[i]<<" with denom norm "<<alphas[i].den().norm()<<endl;
      cout << "sigmas:\n";
      for(int i=0; i<n_sigmas; i++) cout<<i<<": "<<sigmas[i]<<endl;
    }
#ifdef CHECK_TRIANGLES
  int tok = check_triangles(debug);
  assert(tok);
#endif
#ifdef CHECK_SQUARES
  int sok = check_squares(debug);
  assert(sok);
#endif
#ifdef CHECK_HEXAGONS
  int hok = check_hexagons(debug);
  assert(hok);
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
  alpha_denoms.insert(c);
  M_alphas.push_back(M);
  alpha_ind[alpha.coords()] = n_alphas;
  alpha_ind[reduce_to_rectangle(alpha).coords()] = n_alphas;
  n_alphas++;
}

// If r1*r2 = -1 mod s with r1, r2 distinct we have r1/s, -r1/s, r2/s,
// -r2/s with matrices [r2,t;s,-r1] and similar.  When r1==r2 we have
// a "minus pair" with r1^2=-1 (mod s), and when r1==-r2 we have a
// "plus pair" with r1^2=+1 (mod s).

void add_alpha_orbit(const Quad& s, const Quad& r1, const Quad& r2)
{
  Quad t = -(r1*r2+Quad::one);
  assert(div(s,t));
  t /= s;
  if (r1==r2) // "-" pair, r1^2=-1 (mod s)
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
  if (r1==-r2) // "+" pair, r1^2=+1 (mod s)
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
  RatQuad sigma(r,s), msigma(-r,s);
  sigmas.push_back(sigma);
  if (!sigma.is_infinity())
    {
      sigma_ind[sigma.coords()] = n_sigmas;
      sigma_ind[reduce_to_rectangle(sigma).coords()] = n_sigmas;
    }
  if (s.is_zero() || s==TWO) // don't also include -sigma
    {
      sigma_flip.push_back(n_sigmas);     // identity
      n_sigmas+=1;
    }
  else
    {
      sigmas.push_back(msigma);
      sigma_ind[msigma.coords()] = n_sigmas+1;
      sigma_ind[reduce_to_rectangle(msigma).coords()] = n_sigmas+1;
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
  INT D(d);
  Quad w = Quad::w;
  Quad two(2);

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
      Quad u( (D-1)/2 );
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
      Quad u((D-3)/8);  // = 2, 5, 8, 20 for d=19,43,67,163 = 3 (mod 8) so 2 is inert
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
  Quad w = Quad::w;
  Quad three(3);

  switch (d%12) {
  case 1: case 10:
    {
      add_alpha_orbit(three, w, w);  // - pair
      add_alpha_orbit(three, Quad::one+w, Quad::one-w); // 4-some
      // no sigmas (3 inert)
      break;
    }
  case 7:
    {
      if (d>31)
        add_alpha_orbit(three, w, Quad::one-w); // 4-some
      if (d>19) // e.g. 43, 67, 163
        add_alpha_orbit(three, Quad::one+w, Quad::one+w); // - pair
      // no sigmas (3 inert)
      break;
    }
  case 2: case 5:
    {
      if (d!=5)
        {
          add_alpha_orbit(three, w, -w); // + pair
          add_sigma_orbit(Quad::one+w,three);
          add_sigma_orbit(Quad::one-w,three);
        }
      break;
    }
  case 11:
    {
      if (d>=35)
        {
          add_alpha_orbit(three, Quad::one+w, -Quad::one-w); // + pair
          add_sigma_orbit(w,three);
          add_sigma_orbit(w-Quad::one,three);
        }
      break;
    }
  case 3:
    {
      if (d>15)
        {
          add_alpha_orbit(three, w, w-Quad::one); // 4-some
          add_sigma_orbit(Quad::one+w,three);
        }
      break;
    }
  case 6:  case 9:
    {
      if (d>6)
        {
          add_alpha_orbit(three, w+Quad::one, w-Quad::one); // 4-some
          add_sigma_orbit(w,three);
          break;
        }
    }
  } // d%12
}


/****************************************************************

 Code defining triangle relations

***************************************************************/

int check_aaa_triangle(const POLYGON& T, int verbose)
{
  const vector<int>& tri = T.indices;
  const Quad& u = T.shifts[0];
  if (verbose)
    cout<<"Checking aaa-triangle ("<<tri<<","<<u<<")"<<endl;
  mat22 Mi=M_alphas[tri[0]], Mj=M_alphas[tri[1]], Mk=M_alphas[tri[2]];
  RatQuad x = (Mi(Mj.preimage_oo()+u) - Mk.preimage_oo());
  return (x.is_integral());
}

int check_aas_triangle(const POLYGON& T, int verbose)
{
  const vector<int>& tri = T.indices;
  const Quad& u = T.shifts[0];
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

int check_triangles(int verbose)
{
  int all_ok = 1;
  for ( const auto& T : aaa_triangles)
    {
      int ok = check_aaa_triangle(T, verbose);
      assert(ok);
      all_ok = all_ok && ok;
    }
  for ( const auto& T : aas_triangles)
    {
      int ok = check_aas_triangle(T, verbose);
      assert(ok);
      all_ok = all_ok && ok;
    }
  return all_ok;
}

/****************************************************************

 Code defining square relations

***************************************************************/

int check_square(const POLYGON& squ, int verbose)
{
  // Check:  the square {{i,j,k,l},{x,y,z}} has vertices {alpha_i, oo, alpha[j']+z, beta}
  // where beta = z + M_j(x+alpha[k']) = M_i'(y+alpha_l),
  // so that M_i(T^z(M_j(x+alpha[k']))) = y+alpha_l.

  // Edges:

  // {alpha_i, oo} = (I)_i
  // {oo, alpha_j'+z} = (T^z*M_j)_j
  // {alpha_j'+z, beta} = (T^z*M_j*T^x*M_k)_k
  // {beta, alpha_i} = (M_i'*T^y)_l

  const vector<int>& ijkl = squ.indices;  // int i=ijkl[0], j=ijkl[1], k=ijkl[2], l=ijkl[3];
  const vector<Quad>& xyz = squ.shifts;
  if (verbose)
    cout<<"Checking square "<<ijkl<<", "<<xyz<<endl;
  Quad x = xyz[0], y=xyz[1], z=xyz[2];
  mat22 Mi=M_alphas[ijkl[0]], Mj=M_alphas[ijkl[1]], Mk=M_alphas[ijkl[2]], Ml=M_alphas[ijkl[3]];
  RatQuad alpha1 = x + RatQuad(Mk.entry(0,0),Mk.entry(1,0));  // = x+alpha_k'
  RatQuad alpha2 = y + RatQuad(-Ml.entry(1,1),Ml.entry(1,0)); // = y+alpha_l
  mat22 M = Mi*mat22::Tmat(z)*Mj;
  return ((M.entry(0,0)*alpha1+M.entry(0,1))/(M.entry(1,0)*alpha1+M.entry(1,1)) == alpha2);
}

int check_squares(int verbose)
{
  int all_ok = 1;
  for ( const auto& Q : squares)
    {
      int ok = check_square(Q, verbose);
      all_ok = ok && all_ok;
      assert(ok);
    }
  return all_ok;
}

/****************************************************************

 Code defining hexagon relations

***************************************************************/

int check_hexagon(const POLYGON& hex, int verbose)
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

  const vector<int>& ijklmn = hex.indices;
  const vector<Quad>& ux1y1x2y2 = hex.shifts;
  if (verbose)
    cout<<"Checking hexagon "<<ijklmn<<", "<<ux1y1x2y2<<endl;
  Quad u = ux1y1x2y2[0], x1 = ux1y1x2y2[1], y1 = ux1y1x2y2[2], x2 = ux1y1x2y2[3], y2 = ux1y1x2y2[4];
  int i=ijklmn[0], j=ijklmn[1], k=ijklmn[2], l=ijklmn[3], m=ijklmn[4], n=ijklmn[5];
  RatQuad gamma1 = (M_alphas[alpha_inv[i]]*mat22::Tmat(x1)*M_alphas[alpha_inv[k]])(y1+alphas[m]);
  RatQuad gamma2 = (mat22::Tmat(u)*M_alphas[alpha_inv[j]]*mat22::Tmat(x2)*M_alphas[alpha_inv[l]])(y2+alphas[n]);
  return gamma1==gamma2;
}

int check_hexagons(int verbose)
{
  int all_ok = 1;
  for ( const auto& H : hexagons)
    {
      int ok = check_hexagon(H, verbose);
      all_ok = ok && all_ok;
      assert(ok);
    }
  return all_ok;
}

// Return i such that alphas[i]=a, else -1
int alpha_index(const RatQuad& a)
{
  auto s = alpha_ind.find(a.coords());
  if (s==alpha_ind.end())
    return -1;
  return s->second;
}

// Return i and set t such that alphas[i]+t=a, else -1
int alpha_index_with_translation(const RatQuad& a, Quad& t)
{
  auto s = alpha_ind.find(a.coords());
  if (s==alpha_ind.end())
    s = alpha_ind.find(reduce_to_rectangle(a,t).coords());
  if (s==alpha_ind.end())
    return -1;
  int i = s->second;
  (a-alphas[i]).is_integral(t);
  assert (t+alphas[i]==a);
  return i;
}

// Return i such that sigmas[i]=a, else -1
int sigma_index(const RatQuad& a)
{
  auto s = sigma_ind.find(a.coords());
  if (s==sigma_ind.end())
    return -1;
  return s->second;
}

// Return i and set t such that sigmas[i]+t=a, else -1
int sigma_index_with_translation_new(const RatQuad& z, Quad& shift);
int sigma_index_with_translation_old(const Quad& a, const Quad& b, Quad& shift);

int sigma_index_with_translation(const RatQuad& z, Quad& shift)
{
  return sigma_index_with_translation_old(z.num(), z.den(), shift);
}

// Return i and set t such that sigmas[i]+t=a/b, else -1
int sigma_index_with_translation(const Quad& a, const Quad& b, Quad& shift)
{
  return sigma_index_with_translation_old(a, b, shift);
}

int sigma_index_with_translation_old(const Quad& a, const Quad& b, Quad& shift)
{
  int t = 0;
  Quad r, s;
  for ( const auto& sigma : sigmas )
    {
      if (sigma.is_infinity()) // ignore sigma=1/0
        {
          t++;
          continue;
        }
      r=sigma.num(), s=sigma.den(); // sigma = r/s
      shift = mms(a,s,b,r); //a*s-b*r
      if (div(b*s,shift,shift))  // success! NB shift is divided by b*s before returning
        return t;
      t++;
    } // end of loop over singular points sigma
  return -1;
}

int sigma_index_with_translation_new(const RatQuad& z, Quad& shift)
{
  // cout << "Looking up sigma = " << z << endl;
  auto s = sigma_ind.find(z.coords());
  if (s==sigma_ind.end())
    s = sigma_ind.find(reduce_to_rectangle(z,shift).coords());
  if (s==sigma_ind.end())
    return -1;
  int i = s->second;
  (z-sigmas[i]).is_integral(shift);
  assert (shift+sigmas[i]==z);
  // cout << " sigma is sigmas["<<i<<"] + "<<shift<< endl;
  return i;
}

void read_data(int verbose)
{
  ifstream geodata;
  string line;
  int file_d;
  char G;
  POLYGON poly;

  stringstream ss;
  ss << "geodata/geodata_" << Quad::d << ".dat";
  geodata.open(ss.str().c_str());
  if (!geodata.is_open())
    {
      geodata.open("geodata.dat");
      if (verbose)
        cout << "reading from geodata.dat" <<endl;
      if (!geodata.is_open())
        {
          cout << "No geodata file!" <<endl;
          exit(1);
        }
    }
  else
    if (verbose)
      cout << "reading from " << ss.str() <<endl;

  getline(geodata, line);
  while (!geodata.eof())
    {
      parse_geodata_line(line, file_d, G, poly, verbose);
      istringstream input_line(line);
      switch(G) {
      case 'X':
        {
          if (verbose>1)
            cout<<"Skipping rest of file"<<endl;
          break;
        }
      case '%':
        {
          if (verbose>1)
            cout<<"Skipping line: "<<line<<endl;
          break;
        }
      case 'A': // alpha orbit
        {
          add_alpha_orbit(poly.shifts[0], poly.shifts[1], poly.shifts[2]);
          break;
        }
      case 'S': // sigma orbit
        {
          add_sigma_orbit(poly.shifts[0], poly.shifts[1]);
          break;
        }
      case 'T': // aaa-triangle
        {
          aaa_triangles.push_back(poly);
          break;
        }
      case 'U': // aas-triangle
        {
          aas_triangles.push_back(poly);
          break;
        }
      case 'Q': // square
        {
          squares.push_back(poly);
          break;
        }
      case 'H': // hexagon
        {
          hexagons.push_back(poly);
          break;
        }
      default:
        break;
      }
      if (G=='X')
        break; // don't read any more lines
      else
        getline(geodata, line);
    }
  geodata.close();
}

#endif // obsolete code
