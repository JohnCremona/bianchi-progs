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

// indices i,j,k such that M_i(alpha_j)=alpha_k +translation, giving triangle relations
extern vector<vector<int> > triangles;

// [[i,j,k],u] such that M_i(sigma_j+u)=sigma_k +translation, giving triangle relations
extern vector<pair<vector<int>, Quad>> aas_triangles;

// indices i,j,k,l and x,y such that M_j(x+alpha_k') =  M_i'(y+alpha_l), defining a square relation
extern vector<pair<vector<int>, vector<Quad>> > squares;

// indices i,j,k,l,m,n and x1,y1,x2,y2 defining a hexagon relation
extern vector<pair<vector<int>, vector<Quad>> > hexagons;

// initialization function, to be called before constructing face_relations class (for non-Euclidean fields)
void make_faces();

#endif
