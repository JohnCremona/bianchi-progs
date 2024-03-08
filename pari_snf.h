// FILE PARI_SNF.H: declaration of functions for computing invariants of an integer matrix

#if     !defined(_PARI_SNF_H)
#define _PARI_SNF_H      1       //flags that this file has been included

#include <eclib/templates.h>

// Return a list of the Smith Normal Form invariants of a matrix
// (encoded as vector<vector<int>>)
vector<int> invariants(const vector<vector<int>>& M);

#endif
