// FILE PARI_SNF.H: declaration of functions for computing invariants of an integer matrix

#if     !defined(_PARI_SNF_H)
#define _PARI_SNF_H      1       //flags that this file has been included

#include <eclib/templates.h>
#include "int.h"

// Return a list of the Smith Normal Form invariants of a matrix
// (encoded as vector<vector<int>>)
vector<INT> invariants(const vector<vector<int>>& M);

// Return list of entries of the HNF of M (exncoded as above),
// assuming M has size 2 (i.e. represnts a 2xn matrix):
vector<INT> hnf_invariants(const vector<vector<int>>& M);

#endif
