// FILE SWAN_HOM.H: declaration of functions for computing the integral 1-homology

#if     !defined(_SWAN_HOM_H)
#define _SWAN_HOM_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

// Given alphas (and pluspairs, minuspairs, fours), sigmas, faces,
// return the invariants of H_1 as a Z-module in either the SL2 or GL2
// cases or both.

// group=1 for GL2 only, 2 for SL2 only, 3 for both

vector<vector<int>> integral_homology(const vector<CuspList>& faces,
                                      const CuspList& alphas, const CuspList& sigmas,
                                      const vector<vector<Quad>>& pluspairs,
                                      const vector<vector<Quad>>& minuspairs,
                                      const vector<vector<Quad>>& fours,
                                      int group, int debug=0);

// Return the index of an edge {a,b} in the range 0..#alphas+#sigmas-2
int edge_index(const EDGE& e, const CuspList& alphas, const CuspList& sigmas);

// Return the edge boundary matrix M10 (matrix of delta: 1-chains -> 0-chains)
vector<vector<int>> edge_boundary_matrix(const CuspList& alphas, const CuspList& sigmas);

// Return the image under delta of the face, as a vector of length #alphas+#sigmas-1
vector<int> face_boundary_vector(const CuspList& face, const CuspList& alphas, const CuspList& sigmas);

vector<vector<int>> edge_pairings(const vector<CuspList>& faces,
                                  const CuspList& alphas, const CuspList& sigmas,
                                  const vector<vector<Quad>>& pluspairs,
                                  const vector<vector<Quad>>& minuspairs,
                                  const vector<vector<Quad>>& fours,
                                  int GL2);

vector<vector<int>> face_boundaries(const vector<CuspList>& faces,
                                    const CuspList& alphas, const CuspList& sigmas,
                                    int GL2, int debug=0);
// Return the face boundary matrix M21 (matrix of delta: 2-chains -> 1-chains)
vector<vector<int>> face_boundary_matrix(const vector<CuspList>& faces,
                                         const CuspList& alphas, const CuspList& sigmas,
                                         const vector<vector<Quad>>& pluspairs,
                                         const vector<vector<Quad>>& minuspairs,
                                         const vector<vector<Quad>>& fours,
                                         int GL2);

// Given integer matrices (encoded as vector<vector<int>>) of the boundary maps
// M10: 1-chains -> 0-chains (as from edge_boundary_matrix())
// M21: 2-chains -> 1-chains (as from face_boundary_matrix())
// return the invariants of the integral 1-homology

// NB Both matrices are formed by rows, and act on row-vectors on the right
vector<int> homology_invariants(const vector<vector<int>>& M10, const vector<vector<int>>& M21, int debug=0);

void show_invariants(const vector<int>& v);
ostream& operator<<(ostream& os, const vector<vector<int>>& M);

#endif
