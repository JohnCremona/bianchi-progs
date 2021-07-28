// FILE GEOMETRY.H: declaration of functions and associated data for hyperbolic tessation

#if     !defined(_GEOMETRY_H)
#define _GEOMETRY_H      1       //flags that this file has been included

#include <iostream>

#include "mat22.h"

extern int n_alphas;            // Number of alphas.
extern vector<mat22> M_alphas;  // List of matrices M_a  with det(M_a)=1 such that M_a(a)=oo.
extern vector<int> alpha_inv;   // permutation of order 2 swapping a to a' where M_a(oo)=a'
extern vector<int> edge_pairs;  // indices of first of a pair (r/s, -r/s) with r^2=-1 (mod s)
extern vector<int> edge_pairs_plus; // indices of first of a pair (r/s, -r/s) with r^2=+1 (mod s)
extern vector<int> edge_fours;  // indices of first of a 4-tuple (r1,-r1,r2,-r2) of alphas with r1*r2=-1 (mod s)

// indices of alpha such that M_alpha has order 3, giving cyclic triangle relations
extern vector<int> cyclic_triangles;

// indices i,j,k such that M_i(alpha_j)=alpha_k +translation, giving triangle relations
extern vector<vector<int> > triangles;

// indices i,j,k,l and x,y such that M_j(x+alpha_k') =  M_i'(y+alpha_l), defining a square relation
extern vector<pair<vector<int>, vector<Quad>> > squares;

// Singular points: the 0'th is oo so the useful ones are indexed from 1.

extern vector<RatQuad> singular_points;

#endif
