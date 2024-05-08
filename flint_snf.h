#include "int.h"
#include <flint/fmpz_mat.h>

// Given a matrix as vector<vector<int>>, construct a FLINT fmpz_mat

// NB The caller must have called call fmpz_mat_init() on the supplied
// fmpz_mat_t, and must call fmpz_mat_clear() when finished with it.

void make_mat( fmpz_mat_t A, const vector<vector<int>>& M);
void make_mat( fmpz_mat_t A, const vector<vector<INT>>& M);

// Inversely, given a FLINT fmpz_mat, construct a matrix as vector<vector<int>>

void unmake_mat( fmpz_mat_t A, vector<vector<int>>& M);
void unmake_mat( fmpz_mat_t A, vector<vector<INT>>& M);

// compute the rank via FLINT

long rank(const vector<vector<int>>& M);
long rank(const vector<vector<INT>>& M);

// compute the HNF via FLINT

vector<vector<int>> HNF(const vector<vector<int>>& M);
vector<vector<INT>> HNF(const vector<vector<INT>>& M);

// Return a list of the pivotal columns of the HNF of a matrix
// (encoded as vector<vector<int>>) for which the pivots are =1
vector<int> HNF_pivots(const vector<vector<int>>& M);

// compute the SNF via FLINT

void SNF(fmpz_mat_t& S, fmpz_mat_t& A);

