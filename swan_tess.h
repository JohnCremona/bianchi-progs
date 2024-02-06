// FILE SWAN_TESS.H: declaration of functions for finding the tessellation from alphas and sigmas

#if     !defined(_SWAN_TESS_H)
#define _SWAN_TESS_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

// given vertices and edges, fill in faces:
void fill_faces(POLYHEDRON& P, int verbose=0);

// Return all singular polyhedra
vector<POLYHEDRON>
singular_polyhedra(const CuspList& sigmas, const CuspList& alphas, int verbose);

// return a polyhedron (as a list of EDGEs), from the j'th corner Plist[j]
POLYHEDRON
principal_polyhedron(int j, const CuspList& alphas, const H3pointList& Plist,
                     vector<int>& flags, int verbose=0);

// return a list of all principal polyhedra
vector<POLYHEDRON>
principal_polyhedra(const CuspList& alphas, int verbose=0);

// return a list of all polyhedra, principal and singular
vector<POLYHEDRON>
all_polyhedra(const CuspList& alphas, const CuspList& sigmas, int verbose);

// Given a polygon with either all vertices principal or just one
// singular; if there's a singular vertex, rotate to put it at the
// end; then apply M in SL(2,OK) taking the first two vertices to oo
// and a in alist (reduced fundamental alphas).  Set sing to 1 iff singular.
CuspList normalise_polygon( const CuspList& face, const CuspList& alist, const CuspList& sigmas, int& sing);

// Given a polygon, move the first n (default 1) vertices to the end:
CuspList rotate_polygon( const CuspList& face, int n=1);

// Given a polygon, reverse the order of its vertices:
CuspList reverse_polygon( const CuspList& face);

// extract all oriented faces up to SL2-equivalence, returning lists of
// (0) principal triangles {oo, a1, a2} with a1 reduced fundamental
// (1) principal squares   {oo, a1, a2, a3} with a1 reduced fundamental
// (2) principal hexagons  {oo, a1, a2, a3, a4, a5, a6} with a1 reduced fundamental
// (3) singular triangles  {oo, a, s} with a reduced fundamental, s singular
vector<vector<CuspList>> get_faces( const vector<POLYHEDRON>& all_polys,
                                    const CuspList& alphas, const CuspList& sigmas,
                                    int verbose=0);


#endif
