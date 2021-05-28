// FILE EUCLID.H: pseudo-euclidean algorithm and associated data

#if     !defined(_EUCLID_H)
#define _EUCLID_H      1       //flags that this file has been included

#include <iostream>

#include <eclib/arith.h>
#include "quads.h"

extern int n_alphas;            // Number of alphas.
extern vector<mat22> M_alphas;  // List of matrices M_a  with det(M_a)=1 such that M_a(a)=oo.
extern vector<int> alpha_pairs; // permutation of order 2 swapping a to a' where M_a(oo)=a'
void define_alphas();           // Populate M_alphas.

// pseudo-Euclidean step: applies a translation and M_alpha inversion
// to a/b (or column vector [a;b]) reducing b, also multiplying row
// vector [c.d] my M_alpha on the right.  In the Euclidean case, the
// shift is -q where q=a/b (rounded) and the inversion is via
// S=[0,-1;1,0].

// a,b,c,d are changed in place, and on return, t holds the "type"
// (index of alpha which worked)

static int default_t=-1;
void pseudo_euclidean_step(Quad& a, Quad& b, int& t=default_t, Quad& c1=Quad::zero, Quad& d1=Quad::zero, Quad& c2=Quad::zero, Quad& d2=Quad::zero);

// extern mat22 mat22::identity;
// extern mat22 mat22::J;
// extern mat22 mat22::S;
// extern mat22 mat22::T;
// extern mat22 mat22::U;
// extern mat22 mat22::TS;
// extern mat22 mat22::TiS;
// extern mat22 mat22::R;

#endif
