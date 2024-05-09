// FILE SWAN_SIGMAS.H: declaration of singular point functions for Swan's algorithm

#if     !defined(_SWAN_SIGMAS_H)
#define _SWAN_SIGMAS_H      1       //flags that this file has been included

#include <iostream>
#include <set>

#include "swan_utils.h"

CuspList denom_2_sigmas();
CuspList denom_3_sigmas();

// Given an ideal I, return a list of singular points of class [I]
// (one representative for each orbit under integral translations).
CuspList singular_points_in_class(Qideal I, int verbose=0);

// Return a list of lists of singular points in each ideal class.
vector<CuspList> singular_points_by_class();

// Return one list of all singular points.
CuspList singular_points();

// Return sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)
CuspList sort_singular_points(const CuspList& S, int verbose=0);

// Output sorted list of singular points (oo, denom 2, denom 3, larger denoms in +/- pairs)
void output_singular_points(const CuspList& S, int to_file=1, int to_screen=0);

vector<RatQuad> test_singular_points(int output_level=0);

#endif
