// FILE GEOMETRY.H: declaration of functions and associated data for hyperbolic tessation

#if     !defined(_GEOMETRY_H)
#define _GEOMETRY_H      1       //flags that this file has been included

#include <iostream>

#include "mat22.h"

inline ostream& operator<<(ostream& os, const std::set<Quad>& v)
{
  os <<"{ ";
  copy(v.begin(),v.end(), ostream_iterator<Quad>(os, " "));
  os << "}";
  return os;
}

// data for face relations in homology, not initialized by defult

struct POLYGON {
  vector<int> indices; // length n for an n-gon
  vector<Quad> shifts; // length n for an n-gon
};

inline int operator==(const POLYGON& P1, const POLYGON& P2) {
  return P1.indices==P2.indices && P1.shifts==P2.shifts;
}

// reads the polygons (T,U,Q,H) from geodata file: returns 4 lists, of
// aaa-triangles, aas-triangles, squares, hexagons
vector<vector<POLYGON>> read_polygons(string subdir="", int verbose=0);
void parse_geodata_line(const string& line, int& file_d, char& G, POLYGON& poly, int verbose=0);

/***************************** obsolete code using globals ***************/
#if (0)

extern CuspList sigmas; // Singular points: the 0'th is oo, the rest are indexed from 1.
extern int n_sigmas;           // Number of sigmas
extern map<vector<RAT>, int> sigma_ind; // Index of a sigma in the list (from coords as key)

// Base points for principal adges {alpha,oo}
extern CuspList alphas;
extern int n_alphas;                // Number of alphas
extern map<vector<RAT>, int> alpha_ind; // Index of an alpha in the list (from coords as key)
extern std::set<Quad> alpha_denoms; // Denominators of alphas
extern vector<mat22> M_alphas;   // List of M_a with det(M_a)=1 such that M_a(a)=oo and M_a(oo) in alphas

extern vector<int> alpha_inv;    // permutation of order 2 swapping a to a' where M_a(oo)=a'
extern vector<int> alpha_flip;   // permutation of order 2 swapping alpha to -alpha mod 1
extern vector<int> sigma_flip;   // permutation of order 2 swapping sigma to -sigma mod 1

extern vector<int> edge_pairs_minus; // indices of first of a pair (r/s, -r/s) with r^2=-1 (mod s)
extern vector<int> edge_pairs_plus;  // indices of first of a pair (r/s, -r/s) with r^2=+1 (mod s)
extern vector<int> edge_fours;  // indices of first of a 4-tuple (r1,-r1,r2,-r2) of alphas with r1*r2=-1 (mod s)

aaa triangles: [[i,j,k],[u]] such that M_i(alpha_j+u)=alpha_k +translation
aas triangles: [[i,j,k],[u]] such that M_i(sigma_j+u)=sigma_k +translation
squares:  [[i,j,k,l],[x,y,z]] such that M_j(x+alpha_k') =  z + M_i'(y+alpha_l)
hexagons: [[i,j,k,l,m,n], [u,x1,y1,x2,y2]]

extern vector<POLYGON> aaa_triangles;
extern vector<POLYGON> aas_triangles;
extern vector<POLYGON> squares;
extern vector<POLYGON> hexagons;

// Lookup functions for alphas and sigmas

// Return i such that alphas[i]=a, else -1
int alpha_index(const RatQuad& a);
// Return i and set t such that alphas[i]+t=a, else -1
int alpha_index_with_translation(const RatQuad& a, Quad& t);

// Return i such that sigmas[i]=z, else -1
int sigma_index(const RatQuad& z);
// Return i and set t such that sigmas[i]+t=z, else -1
int sigma_index_with_translation(const RatQuad& z, Quad& t);
// Return i and set t such that sigmas[i]+t=a/b, else -1
int sigma_index_with_translation(const Quad& a, const Quad& b, Quad& t);

// The following require M_alphas to be defined properly, and for aas triangles also that sigmas is defined.
int check_aaa_triangle(const POLYGON& T, int verbose=0);
int check_aas_triangle(const POLYGON& T, int verbose=0);
int check_square(const POLYGON& squ, int verbose=0);
int check_hexagon(const POLYGON& hex, int verbose=0);

#endif // obsolete code

#endif
