// FILE GEOMETRY.H: declaration of functions and associated data for hyperbolic tessation

#if     !defined(_GEOMETRY_H)
#define _GEOMETRY_H      1       //flags that this file has been included

#include <iostream>

#include "mat22.h"

extern vector<RatQuad> sigmas; // Singular points: the 0'th is oo, the rest are indexed from 1.
extern int n_sigmas;           // Number of sigmas
extern map<vector<RAT>, int> sigma_ind; // Index of a sigma in the list (from coords as key)
extern std::set<Quad, Quad_comparison> sigma_denoms; // Denominators of finite singular points

// Base points for principal adges {alpha,oo}
extern vector<RatQuad> alphas;
extern int n_alphas;                // Number of alphas
extern map<vector<RAT>, int> alpha_ind; // Index of an alpha in the list (from coords as key)
extern std::set<Quad, Quad_comparison> alpha_denoms; // Denominators of alphas
extern vector<mat22> M_alphas;   // List of M_a with det(M_a)=1 such that M_a(a)=oo and M_a(oo) in alphas

extern vector<int> alpha_inv;    // permutation of order 2 swapping a to a' where M_a(oo)=a'
extern vector<int> alpha_flip;   // permutation of order 2 swapping alpha to -alpha mod 1
extern vector<int> sigma_flip;   // permutation of order 2 swapping sigma to -sigma mod 1

extern vector<int> edge_pairs_minus; // indices of first of a pair (r/s, -r/s) with r^2=-1 (mod s)
extern vector<int> edge_pairs_plus;  // indices of first of a pair (r/s, -r/s) with r^2=+1 (mod s)
extern vector<int> edge_fours;  // indices of first of a 4-tuple (r1,-r1,r2,-r2) of alphas with r1*r2=-1 (mod s)

inline ostream& operator<<(ostream& os, const std::set<Quad, Quad_comparison>& v)
{
  os <<"{ ";
  copy(v.begin(),v.end(), ostream_iterator<Quad>(os, " "));
  os << "}";
  return os;
}

// data for face relations in homology, not initialized by defult

// indices of alpha such that M_alpha has order 3, giving cyclic triangle relations
extern vector<int> cyclic_triangles;

// typedefs for readability.

typedef pair<vector<int>, Quad> TRIANGLE;
typedef pair<vector<int>, vector<Quad>> POLYGON;

// aaa triangles: [[i,j,k],u] such that M_i(alpha_j+u)=alpha_k +translation
// aas triangles: [[i,j,k],u] such that M_i(sigma_j+u)=sigma_k +translation
// (We could make TRIANGLE a special case of POLYGON where the second component is a vector of length 1.)

// squares:  [[i,j,k,l],[x,y,z]] such that M_j(x+alpha_k') =  z + M_i'(y+alpha_l)
// hexagons: [[i,j,k,l,m,n], [u,x1,y1,x2,y2]]

extern vector<TRIANGLE> aaa_triangles;
extern vector<TRIANGLE> aas_triangles;
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

#endif
