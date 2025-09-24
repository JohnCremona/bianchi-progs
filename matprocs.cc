// FILE MATPROCS.CC:  implementation of interface with NTL mat_ZZ and related classes

#include <eclib/linalg.h>
#include "matprocs.h"

factor_comparison fact_cmp;
poly_comparison poly_cmp;

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

// return f(X/c)*c^d: multiply coeff(f,i) by c^(d-i)
ZZX scale_poly_up(const ZZX& f, const ZZ& c)
{
  if (c==1) return f;
  ZZX g = f;
  ZZ cpow(1);
  long d = deg(f);
  for (int i=0; i<=d; i++)
    {
      SetCoeff(g, d-i, cpow*coeff(g, d-i));
      if (i<d)
        cpow *= c;
    }
  return g;
}

// return f(c*X)/c^d: divide coeff(f,i) by c^(d-i)
// NB only use when divisions are exact
ZZX scale_poly_down(const ZZX& f, const ZZ& c)
{
  if (c==1) return f;
  ZZX g = f;
  ZZ cpow(1);
  long d = deg(f);
  for (int i=0; i<=d; i++)
    {
      SetCoeff(g, d-i, coeff(g, d-i)/cpow);
      if (i<d)
        cpow *= c;
    }
  return g;
}

// return f(X) mod m
ZZX reduce_poly(const ZZX& f, const ZZ& m)
{
  if (m==0) return f;
  ZZX g = f;
  long d = deg(f);
  for (int i=0; i<=d; i++)
    SetCoeff(g, d-i, mod(coeff(g, d-i), m));
  return g;
}

// Squarefree test
int IsSquareFree(const ZZX& f)
{
  return deg(GCD(f, diff(f))) == 0;
}

ZZX scaled_charpoly(const mat_ZZ& A, const ZZ& den, const scalar& m)
{
  ZZX charpol;
  CharPoly(charpol, A);
  return reduce_poly(scale_poly_down(charpol, den), to_ZZ(m));
}

// return A mod m (or just A if m==0)
mat_ZZ reduce_mat(const mat_ZZ& A, const ZZ& m)
{
  if (m==0) return A;
  int nr= A.NumRows(), nc = A.NumCols();
  mat_ZZ B = A;
  for (int i=1; i<=nr; i++)
    for (int j=1; j<=nc; j++)
      B(i,j) = mod(B(i,j), m);
  return B;
}

// evaluate f(A) (assumes f monic)
mat_ZZ evaluate(const ZZX& f, const mat_ZZ& A)
{
  long d = deg(f);
  mat_ZZ fA = A, I;
  ident(I, A.NumRows());
  for(int i=d-1; i>=0; i--)
    {
      fA += coeff(f,i)*I;
      if(i)
        fA *= A;
    }
  return fA;
}

// evaluate f(A) (assumes f monic)
mat evaluate(const ZZX& f, const mat& A)
{
  long d = deg(f);
  mat_m mA = to_mat_m(A);
  mat_m fA(mA);
  for(int i=d-1; i>=0; i--)
    {
      fA = addscalar(fA,coeff(f,i));
      if(i)
        fA = fA*mA;
    }
  return to_mat(fA);
}

int check_involution(const mat_ZZ& A, scalar den, const scalar& m, int verbose)
{
  int res = IsDiag(reduce_mat(sqr(A), to_ZZ(m)), A.NumRows(), to_ZZ(den*den));
  if (verbose)
    cout << (res? "Involution!": "NOT an involution....") << "\n";
  return res;
}

int commute(const mat_ZZ& A, const mat_ZZ& B, const scalar& m)
{
  return IsZero(reduce_mat(A*B-B*A, to_ZZ(m)));
}

// check that a matrix commutes with all those in a list:
int check_commute(const mat_ZZ& A, const vector<mat_ZZ>& Blist, const scalar& modulus)
{
  return std::all_of(Blist.begin(), Blist.end(), [A, modulus] (const mat_ZZ& B) {return commute(A,B,modulus);});
}

// display factors of a polynomaial:
void display_factors(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
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
