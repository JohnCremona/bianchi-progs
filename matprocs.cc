// FILE MATPROCS.CC:  implementation of interface with NTL mat_ZZ and related classes

#include <eclib/linalg.h>
#include "matprocs.h"

template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A)
{
  int d = A.nrows();
  // copy into an NTL matrix:
  mat_ZZ ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j)=NTL::conv<ZZ>(A(i,j));
  return ntl_A;
}
template mat_ZZ mat_to_mat_ZZ<int>(Zmat<int> A);
template mat_ZZ mat_to_mat_ZZ<long>(Zmat<long> A);
template mat_ZZ mat_to_mat_ZZ<ZZ>(Zmat<ZZ> A);

mat_ZZ_p mat_to_mat_ZZ_p(mat A)
{
  int d = A.nrows();

  // copy into an NTL matrix:
  mat_ZZ_p ntl_A;
  ntl_A.SetDims(d,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      ntl_A(i,j)=NTL::conv<ZZ_p>(A(i,j));
  return ntl_A;
}

// compute char poly of A:
ZZX charpoly(const mat_ZZ& A)
{
  ZZX charpol;
  CharPoly(charpol, A);
  return charpol;
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
template<class T>
mat_m evaluate(const ZZX& f, const Zmat<T>& A)
{
  long d = deg(f);
  mat_m mA(to_mat_m(A));
  mat_m fA(mA);
  for(int i=d-1; i>=0; i--)
    {
      fA = addscalar(fA, coeff(f,i));
      if(i)
        fA = fA*mA;
    }
  return fA;
}

template mat_m evaluate<int>(const ZZX& f, const Zmat<int>& A);
template mat_m evaluate<long>(const ZZX& f, const Zmat<long>& A);
template mat_m evaluate<ZZ>(const ZZX& f, const Zmat<ZZ>& A);

// p should be monic:
mat_ZZ CompanionMatrix(const ZZX& p)
{
  int d = deg(p);
  mat_ZZ A;
  A.SetDims(d,d);
  ZZ one(1);
  for(int i=1; i<d; i++)
    {
      A(i+1,i) = one;
      A(i,d) = -coeff(p, i-1);
    }
  A(d,d) = -coeff(p, d-1);
  return A;
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

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vec_m& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i+1];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

// Linear combinarion of n>0 matrices, all dxd
mat_m lin_comb_mats(const vector<ZZ>& co, const vector<mat_m>& mats)
{
  int n = mats.size(), d = mats[0].nrows();
  mat_m a(d,d);
  for (int i=0; i<n; i++)
    {
      ZZ c = co[i];
      if (c!=0)
        a += c*mats[i];
    }
  return a;
}

// same as m.output(cout) except no newlines between rows
template<class T>
void output_flat_matrix(const Zmat<T>& m, ostream&s)
{
  vector<T> entries = m.get_entries();
  auto mij = entries.begin();
  s << "[";
  long nr = m.nrows();
  while(nr--)
    {
      long nc = m.ncols();
      s<<"[";
      while(nc--)
        {
          s<<(*mij++);
          if(nc)
            s<<",";
        }
      s<<"]";
      if(nr)
        s<<",";
    }
  s << "]";
}

template void output_flat_matrix<int>(const Zmat<int>& m, ostream&s);
template void output_flat_matrix<long>(const Zmat<long>& m, ostream&s);
template void output_flat_matrix<ZZ>(const Zmat<ZZ>& m, ostream&s);

// rank of an NTL matrix:
long rank(mat_ZZ A)
{
  ZZ d2;
  return NTL::image(d2, A);
}

// nullity of an NTL matrix:
long nullity(mat_ZZ A)
{
  return A.NumRows()-rank(A);
}
