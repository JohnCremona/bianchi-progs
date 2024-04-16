// FILE MATPROCS.CC:  implementation of interface with NTL mat_ZZ and related classes

#include <eclib/interface.h>
#include <eclib/arith.h>
#include <eclib/matrix.h>
#include <eclib/mmatrix.h>
#include "matprocs.h"

factor_comparison fact_cmp;

mat_ZZ mat_to_mat_ZZ(mat A)
{
  int i, j, d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ>(A(i,j));
  return ntl_A;
}

mat_ZZ_p mat_to_mat_ZZ_p(mat A)
{
  int i, j, d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ_p ntl_A;
  ntl_A.SetDims(d,d);
  for(i=1; i<=d; i++)
    for(j=1; j<=d; j++)
      ntl_A(i,j)=conv<ZZ_p>(A(i,j));
  return ntl_A;
}

ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den)
{
  ZZX charpol;
  CharPoly(charpol, A);
  long d = deg(charpol);
  if (den>1)
    {
      bigint dpow(1);
      for(int i=0; i<=d; i++)
        {
          SetCoeff(charpol, d-i, coeff(charpol, d-i)/dpow);
          dpow *= den;
        }
    }
  return charpol;
}

int check_involution(const mat_ZZ& A, long den, int verbose)
{
  int res = IsDiag(power(A,2), A.NumRows(), to_ZZ(den*den));
  if (verbose)
    cout << (res? "Involution!": "NOT an involution....") << "\n";
  return res;
}

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist)
{
  return std::all_of(Blist.begin(), Blist.end(), [A] (const mat_ZZ& B) {return A*B==B+A;});
}

// display factors of a polynomaial:
void display_factors(const ZZX& f)
{
  ZZ content; vec_pair_ZZX_long factors;
  factor(content, factors, f);
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
