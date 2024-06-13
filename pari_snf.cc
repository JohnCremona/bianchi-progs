// FILE PARI_SNF.CC: implementation of functions for computing invariants of an integer matrix

#include "pari_snf.h"
#include <eclib/interface.h> // for getenv_with_default
#include <pari/pari.h>

// This is only the default of the environmant variable PARI_SIZE is not set
#define DEFAULT_PARI_SIZE 10000000000 // 10^10 = 10GB approx (10^9 not enough for d=911)
#define DEFAULT_PARI_MAX_PRIME 1000000

void eclib_pari_init(long max_prime=DEFAULT_PARI_MAX_PRIME)
{
  if (!avma) {
    long pari_size = strtol(getenv_with_default("PARI_SIZE", "DEFAULT_PARI_SIZE").c_str(), NULL, 0);
    if (pari_size==0) // e.g. syntax error in the environment variable PARI_SIZE
      pari_size = DEFAULT_PARI_SIZE;
#ifdef DEBUG_GPFACT
    std::cout<<"calling pari_init with pari_size = "<<pari_size<<endl;
#endif
    // the first parameter is the maximum stack size in bytes
    // the second parameter is the maximum precomputed prime
    pari_init(pari_size, max_prime);
  }
}

vector<int> invariants(const vector<vector<int>>& M)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer

  long nrows = M.size(), ncols = M[0].size();
  // cout << "\nnrows="<<nrows<<", ncols="<<ncols<<endl;
  GEN A = zeromatcopy(ncols, nrows);
  // cout << "created A"<<endl;
  for (int i=0; i<nrows; i++)
    for (int j=0; j<ncols; j++)
      gcoeff(A, j+1, i+1) = stoi(M[i][j]);
  // cout << "filled A"<<endl;
  GEN S = ZM_snf(A);
  // cout << "computed S"<<endl;
  long s = lg(S)-1; // itos(gel(matsize(S), 2));
  // cout << "computed size of S = "<<s<<endl;
  vector<int> invs;
  for (int i=0; i<s; i++)
    {
      int d = itos(gel(S,s-i)); // reversing order
      if (d!=1)
        invs.push_back(d);
    }
  avma=av;
  return invs;
}

vector<long> hnf_invariants(const vector<vector<int>>& M)
{
  // assuming nrows = 2
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer

  long nrows = M.size(), ncols = M[0].size();
  // cout << "\nnrows="<<nrows<<", ncols="<<ncols<<endl;
  GEN A = zeromatcopy(ncols, nrows);
  // cout << "created A"<<endl;
  for (int i=0; i<nrows; i++)
    for (int j=0; j<ncols; j++)
      gcoeff(A, j+1, i+1) = stoi(M[i][j]);
  // cout << "filled A"<<endl;
  GEN H = ZM_hnf(A);
  // cout << "computed H"<<endl;
  vector<long> invs = {itos(gcoeff(H,1,1)), itos(gcoeff(H,1,2)), itos(gcoeff(H,2,1)), itos(gcoeff(H,2,2))};
  avma=av;
  return invs;
}
