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

// Given a polygon, negate its vertices (i.e. apply a transformation in GL2 not SL2):
CuspList negate_polygon( const CuspList& face);

// extract all oriented faces up to GL2-equivalence, rotation and reflection,
// returning a single mixed list of:
// - principal triangles {a1, oo, a2} with a1 reduced fundamental
// - principal squares   {a1, oo, a2, a3} with a1 reduced fundamental
// - principal hexagons  {a1, oo, a2, a3, a4, a5, a6} with a1 reduced fundamental
// - singular triangles  {a, oo, s} with a reduced fundamental, s singular

// M32 returns a matrix (encoded as vector<vector<int>>) with one row
// per polyhedron giving its boundary as a Z-linear combination of
// oriented faces

vector<CuspList> get_faces( const vector<POLYHEDRON>& all_polys,
                            const CuspList& alphas, const CuspList& sigmas,
                            vector<vector<int>>& M32,
                            int verbose=0);

// Return complete string encoding one face
string encode_int_list(char type, const vector<INT> data);

// Return string for POLYGON representing an aaa-triangle, aas-triangle, quadrilateral or hexagon
string polygon_string(const POLYGON& P, int sing);

// For any face, return the string which encodes it in the geodata files
string face_encode(const CuspList& face, const CuspList& alphas, const CuspList& sigmas);

// Convert an actual polygon (aaa- or aas-triangle, quadrilateral or
// hexagon) as list of vertices to the POLYGON {{i,j,k},{u}}, setting
// sing to 1 for an aas-triangle, else to 0
POLYGON make_polygon(const CuspList& face, const CuspList& alphas, const CuspList& sigmas, int& sing);

// Convert an actual aaa- or aas-triangle as list of vertices [a,oo,b]
// or [a,oo,s] to a POLYGON {{i,j,k},{u}}, setting sing to 1 for an
// aas-triangle, else to 0
POLYGON make_triangle(const CuspList& T, const CuspList& alphas, const CuspList& sigmas, int& sing);

// Convert an actual quadrilateral as list of vertices [a,oo,b,c] to
// the POLYGON {{i,j,k},{x,y,z}}
POLYGON make_quadrilateral(const CuspList& Q, const CuspList& alphas);

// Convert an actual hexagon as list of vertices [a_i, oo, a_j, b_2,
// gamma, b_1] to the POLYGON {{i,j,k,l,m,n},{u,x1,y1,x2,y2}}
POLYGON make_hexagon(const CuspList& H, const CuspList& alphas);

void output_faces( const vector<vector<CuspList>>& aaa_squ_hex_aas,
                   const CuspList& alphas, const CuspList& sigmas,
                   int to_file, int to_screen);

#endif
