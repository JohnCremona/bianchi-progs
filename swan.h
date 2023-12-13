// FILE SWAN.H: declaration of Swan's algorithm functions

#if     !defined(_SWAN_H)
#define _SWAN_H      1       //flags that this file has been included

#include <iostream>

#include "quads.h"

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).

vector<RatQuad> singular_points_in_class(Qideal I, int verbose=0);

// Return a list of lists of singular points in each ideal class.

vector<vector<RatQuad>> singular_points_by_class();

// Return one list of all singular points.

vector<RatQuad> singular_points();

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

vector<RatQuad> sort_singular_points(const vector<RatQuad> S, int verbose=0);

// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)

void output_singular_points(const vector<RatQuad> S, int to_file=1, int to_screen=0);

#endif
