// FILE MATPROCS.CC:  implementation of interface with NTL mat_ZZ and related classes

#include <eclib/linalg.h>
#include "matprocs.h"

factor_comparison fact_cmp;

mat_ZZ mat_to_mat_ZZ(mat A)
{
  int d = A.nrows();
  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ>(A(i,j));
  return ntl_A;
}

mat_ZZ_p mat_to_mat_ZZ_p(mat A)
{
  int d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ_p ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ_p>(A(i,j));
  return ntl_A;
}

ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den, const ZZ& modulus)
{
  ZZX charpol;
  CharPoly(charpol, A);
  long d = deg(charpol);
  if ((den>1)||(modulus!=0))
    {
      bigint dpow(1);
      for(int i=0; i<=d; i++)
        {
          bigint c = coeff(charpol, d-i)/dpow;
          if (modulus!=0)
            {
              c = mod(c, modulus);
            }
          SetCoeff(charpol, d-i, c);
          dpow *= den;
        }
    }
  return charpol;
}

int check_involution(const mat_ZZ& A, scalar den, const ZZ& modulus, int verbose)
{
  int n = A.NumRows();
  mat_ZZ Asq = sqr(A);
  if (modulus!=0)
    {
      for (int i=1; i<=n; i++)
        for (int j=1; j<=n; j++)
          Asq(i,j) = mod(Asq(i,j), modulus);
    }
  int res = IsDiag(Asq, n, to_ZZ(den*den));
  if (verbose)
    cout << (res? "Involution!": "NOT an involution....") << "\n";
  return res;
}

int commute(const mat_ZZ& A, const mat_ZZ& B, const ZZ& modulus)
{
  mat_ZZ C = A*B-B*A;
  int n = C.NumRows();
  for (int i=1; i<=n; i++)
    for (int j=1; j<=n; j++)
      {
        ZZ Cij = C(i,j);
        if (modulus!=0)
          Cij = mod(Cij, modulus);
        if (Cij!=0) return 0;
      }
  return 1;
}

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const ZZ& modulus)
{
  return std::all_of(Blist.begin(), Blist.end(), [A, modulus] (const mat_ZZ& B) {return commute(A,B,modulus);});
}

// display factors of a polynomaial:
void display_factors(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
  cout<<"Factors of characteristic polynomial are:"<<endl;
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout<<(i+1)<<":\t"<<factors[i].a
          <<"\t(degree "<<deg(factors[i].a)<<")";
      cout<<"\t to power "<<factors[i].b;
      cout<<endl;
    }
}

// rank of an NTL matrix:
long rank(mat_ZZ A)
{
  ZZ d2;
  return image(d2, A);
}

// nullity of an NTL matrix:
long nullity(mat_ZZ A)
{
  return A.NumRows()-rank(A);
}
