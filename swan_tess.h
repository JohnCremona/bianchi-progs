// FILE SWAN_TESS.H: declaration of functions for finding the tessellation from alphas and sigmas

#if     !defined(_SWAN_TESS_H)
#define _SWAN_TESS_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

// given vertices and edges, fill in faces:
void fill_faces(POLYHEDRON& P, int verbose=0);

// return a tetrahedron from a list of its vertices
POLYHEDRON tetrahedron(const CuspList& V);

// return a list of tetrahedra (i.e. lists of 4 cusps (oo, sigmas[j],
// a1, a2) with a1,a2 fundamental); as a side-effect set flags[j]=1
// and flags[j']=1 where M_a_i(sigma[j])=sigma[j'] for i=1,2, for
// each.
vector<POLYHEDRON>
singular_tetrahedra(int j, const CuspList& sigmas, const CuspList& alphas, vector<int>& flags, int verbose=0);

// return a list of all singular tetrahedra
vector<POLYHEDRON>
singular_tetrahedra(const CuspList& sigmas, const CuspList& alphas, int verbose=0);

// return a polyhedron (as a list of EDGEs), from the j'th corner Plist[j]
POLYHEDRON
principal_polyhedron(int j, const CuspList& alphas, const H3pointList& Plist,
                     vector<int>& flags, int verbose=0);

// return a list of all principal polyhedra
vector<POLYHEDRON>
principal_polyhedra(const CuspList& alphas, int verbose=0);

#endif
