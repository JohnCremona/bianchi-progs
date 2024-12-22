// FILE PARI_SNF.CC: implementation of functions for computing invariants of an integer matrix

#include "pari_snf.h"
#include <eclib/interface.h> // for getenv_with_default
#include <pari/pari.h>
#include <assert.h>

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

// Convert a pari t_INT to an INT
const long BLOCK_SIZE = 32;
INT convert_t_INT_to_INT(GEN n, int debug=0)
{
  // if n fits in a long int it is easy:
  if (!is_bigint(n))
    return INT(itos(n));

  // if n<0 convert |n| and then negate:
  if (signe(n)<0)
    return -convert_t_INT_to_INT(negi(n));

  // For positive n > 2^63, write in a small enough base so that the
  // digits can be converted to long ints:
  if(debug)
    pari_printf("Converting t_INT %Ps\n", n);
  GEN digits = binary_2k(n, BLOCK_SIZE); // n's digits in base 2^32, as a t_VEC
  if(debug)
    pari_printf("Digits are %Ps\n", digits);
  INT ans(0);
  for (int i=1; i<lg(digits); i++)
    {
      if (i)
        ans <<= BLOCK_SIZE;
      ans += itos(gel(digits,i));
    }
  if(debug)
    cout << "Returning "<<ans<<endl;
  return ans;
}

void test_convert(int debug=0)
{
  if(debug)
    cout << "In test_convert()" << endl;
  pari_sp av=avma;  // store pari stack pointer
  long n = 461018427387914;
  // long n = 4611686018427387914; // 2**62 + 10
  INT N(n);
  if(debug)
    cout << "n = " << n << endl;
  GEN a = stoi(n);
  if(debug)
    pari_printf("Testing conversion of a = %Ps\n", a);
  assert (convert_t_INT_to_INT(a) == N);
  if(debug)
    cout<<N<<" converted ok"<<endl;
  GEN a10 = addsi(10,a);
  if(debug)
    pari_printf("a10 = a+10 = %Ps\n", a10);
  GEN b = mulii(a,a10);
  if(debug)
    pari_printf("b = a(a+10) = %Ps\n", b);
  INT M = N*(10+N);
  if(debug)
    pari_printf("Testing conversion of %Ps\n", b);
  assert (convert_t_INT_to_INT(b) == M);
  if(debug)
    cout<<M<<" converted ok"<<endl;
  GEN b10 = addsi(10,b);
  GEN c = mulii(b,b10);
  M = M*(10+M);
  if(debug)
    pari_printf("Testing conversion of %Ps\n", c);
  assert (convert_t_INT_to_INT(c) == M);
  if(debug)
    cout<<M<<" converted ok"<<endl;
  avma=av;
}

vector<INT> invariants(const vector<vector<int>>& M)
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
  vector<INT> invs;
  // pari_printf("Invariants in libpari: %Ps\n", S);
  for (int i=0; i<s; i++)
    {
      GEN e = gel(S,s-i); // reversing order
      INT d = convert_t_INT_to_INT(e);
      if (!is_one(d))
        {
          // if (d>2)
          //   {
          //     pari_printf("%Ps converts to ", e);
          //     cout<<d<<endl;
          //   }
          invs.push_back(INT(d));
        }
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
      INT d = convert_t_INT_to_INT(e);
      if (!is_one(d))
        {
          // if (d>2)
          //   {
          //     pari_printf("%Ps converts to ", e);
          //     cout<<d<<endl;
          //   }
          invs.push_back(INT(d));
        }
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
    convert_t_INT_to_INT(gcoeff(H,1,1)),
    convert_t_INT_to_INT(gcoeff(H,1,2)),
    convert_t_INT_to_INT(gcoeff(H,2,1)),
    convert_t_INT_to_INT(gcoeff(H,2,2))};
  avma=av;
  return invs;
}
