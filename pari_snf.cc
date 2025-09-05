// FILE PARI_SNF.CC: implementation of functions for computing invariants of an integer matrix

#include "pari_snf.h"
#include <eclib/interface.h> // for getenv_with_default
#include <eclib/pari_init.h>
#include <eclib/convert.h>
#include <assert.h>

using PARI::zeromatcopy;
using PARI::ZM_rank;
using PARI::ZM_inv;
using PARI::ZM_mul;
using PARI::ZM_snf;
using PARI::ZM_hnf;
using PARI::ZM_hnfall;
using PARI::RgM_dimensions;
using PARI::rowslice;

// Convert a pari t_INT to an INT
INT PARI_to_INT(GEN n)
{
  return INT(*PARI_to_FLINT(n));
}

vector<INT> invariants(const vector<vector<int>>& M)
{
  eclib_pari_init();

  pari_sp av=avma;  // store pari stack pointer

  long nrows = M.size(), ncols = M[0].size();
  GEN A = zeromatcopy(ncols, nrows);
  for (int i=0; i<nrows; i++)
    for (int j=0; j<ncols; j++)
      gcoeff(A, j+1, i+1) = stoi(M[i][j]);
  GEN S = ZM_snf(A);
  long s = lg(S)-1; // itos(gel(matsize(S), 2));
  vector<INT> invs;
  for (int i=0; i<s; i++)
    {
      GEN e = gel(S,s-i); // reversing order
      INT d = PARI_to_INT(e);
      if (!is_one(d))
        invs.push_back(INT(d));
    }
  avma=av;
  return invs;
}

vector<INT> homology_invariants_via_pari(const vector<vector<int>>& M10, const vector<vector<int>>& M21, int debug)
{
  eclib_pari_init();
  pari_sp av=avma;  // store pari stack pointer
  long nr, nc;

  // M10 represents a n1xn0 matrix and M21 a n2xn1, with M21*M10=0
  long n0 = M10[0].size(), n1 = M10.size(), n2 = M21.size();
  assert (n1==(long)M21[0].size());

  GEN A10 = zeromatcopy(n0, n1);
  for (int i=0; i<n0; i++)
    for (int j=0; j<n1; j++)
      gcoeff(A10, i+1, j+1) = stoi(M10[j][i]);

  if (debug)
    cout<<"Created transposed pari matrix A10 of size "<<n0<<"x"<<n1<<endl;

  GEN A21 = zeromatcopy(n1, n2);
  for (int i=0; i<n1; i++)
    for (int j=0; j<n2; j++)
      gcoeff(A21, i+1, j+1) = stoi(M21[j][i]);

  if (debug)
    cout<<"Created transposed pari matrix A21 of size "<<n1<<"x"<<n2<<endl;

  GEN U = zeromatcopy(n1, n1);
  GEN H = ZM_hnfall(A10, &U, 0); // remove=0
  long r = ZM_rank(H);
  if (debug)
    cout<<"HNF(A10) -> U, H where rank(H)="<<r<<endl;
  GEN den = stoi(1);
  GEN Uinv = ZM_inv(U, &den);
  if (debug)
    {
      cout<<"Computed U^{-1};";
      pari_printf(" denom = %Ps", den);
      cout<<endl;
    }
  GEN M = ZM_mul(Uinv, A21);
  if (debug)
    {
      RgM_dimensions(M,&nr,&nc);
      cout<<"U^{-1}*A21 -> M ("<<nr<<"x"<<nc<<")"<<endl;
    }
  GEN A = rowslice(M, 1, n1-r);
  if (debug)
    {
      RgM_dimensions(A,&nr,&nc);
      cout << "M with last r rows dropped -> A ("<<nr<<"x"<<nc<<")"<<endl;
      cout << "Computing S = SNF(A)..."<<endl;
    }
  GEN S = ZM_snf(A);
  if (debug) cout << "...done."<<endl;
  long s = lg(S)-1; // itos(gel(matsize(S), 2));
  if (debug) cout << " size of S = "<<s<<endl;
  vector<INT> invs;
  if (debug>1) pari_printf("Invariants in libpari: %Ps\n", S);
  for (int i=0; i<s; i++)
    {
      GEN e = gel(S,s-i); // reversing order
      INT d = PARI_to_INT(e);
      if (!is_one(d))
        invs.push_back(INT(d));
    }
  avma=av;
  if (debug) cout<<"non-trivial invariants:  " << invs <<endl;
  return invs;
}

vector<INT> hnf_invariants(const vector<vector<int>>& M)
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
  vector<INT> invs = {
    PARI_to_INT(gcoeff(H,1,1)),
    PARI_to_INT(gcoeff(H,1,2)),
    PARI_to_INT(gcoeff(H,2,1)),
    PARI_to_INT(gcoeff(H,2,2))};
  avma=av;
  return invs;
}
