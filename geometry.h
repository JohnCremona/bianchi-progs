// FILE GEOMETRY.H: declaration of functions and associated data for hyperbolic tessation

#if     !defined(_GEOMETRY_H)
#define _GEOMETRY_H      1       //flags that this file has been included

#include <iostream>

#include "mat22.h"

extern vector<RatQuad> sigmas; // Singular points: the 0'th is oo, the rest are indexed from 1.
extern int n_sigmas;           // Number of sigmas

// Base points for principal adges {alpha,oo}
extern vector<RatQuad> alphas;
extern int n_alphas;            // Number of alphas
extern vector<mat22> M_alphas;  // List of matrices M_a  with det(M_a)=1 such that M_a(a)=oo.

extern vector<int> alpha_inv;   // permutation of order 2 swapping a to a' where M_a(oo)=a'
extern vector<int> alpha_flip;   // permutation of order 2 swapping alpha to -alpha mod 1
extern vector<int> sigma_flip;   // permutation of order 2 swapping sigma to -sigma mod 1

extern vector<int> edge_pairs_minus; // indices of first of a pair (r/s, -r/s) with r^2=-1 (mod s)
extern vector<int> edge_pairs_plus;  // indices of first of a pair (r/s, -r/s) with r^2=+1 (mod s)
extern vector<int> edge_fours;  // indices of first of a 4-tuple (r1,-r1,r2,-r2) of alphas with r1*r2=-1 (mod s)

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

#endif
