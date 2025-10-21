// FILE MATPROCS.CC:  implementation of interface with NTL mat_ZZ and related classes

#include <eclib/linalg.h>
#include "matprocs.h"

factor_comparison fact_cmp;
poly_comparison poly_cmp;
factor_modp_comparison fact_modp_cmp;
poly_modp_comparison poly_modp_cmp;

// pretty output for integer polynomials
string monomial_string(int i, const string& var)
{
  ostringstream s;
  if (i>0) s << var;
  if (i>1) s << "^" << i;
  return s.str();
}

string polynomial_string(const vector<ZZ>& coeffs, const string& var)
{
  //  cout<<"\npolynomial_string("<<coeffs<<")\n";
  if (std::all_of(coeffs.begin(), coeffs.end(), [](const ZZ&c){return IsZero(c);}))
    return "0";
  int d = coeffs.size()-1;
  ZZ c;
  ostringstream s;
  if (d==0)
    {
      s << coeffs[0];
      return s.str();
    }
  // All non-constant terms:
  for (int i=d; i>0; i--)
    {
      c = coeffs[i];
      if (c==0)
        continue;
      if (c>1)
        {
          if (i<d) // no + needed on leading term
            s << "+";
          s << c << "*";
        }
      else if (c==1 && i<d) s << "+";
      else if (c==-1) s << "-";
      else if (c<-1) s << "-" << abs(c) << "*";
      s << monomial_string(i, var);
      //cout<<" - after i="<<i<<": "<<s.str()<<endl;
    }
  // Constant term:
  c = coeffs[0];
  if (c>0) s << "+" << c;
  else if (c<0) s << "-" <<abs(c);
  //cout<<" - after i=0: "<<s.str()<<endl;

  return s.str();
}

string polynomial_string(const vec_m& coeffs, const string& var)
{
  if (trivial(coeffs))
    return "0";
  int d = dim(coeffs); // one less than 'degree'
  vector<ZZ> co(d);
  for(int i=0; i<d; i++)
    co[i] = coeffs[i+1];
  return polynomial_string(co, var);
}

vector<ZZ> coeffs(const ZZX& p)
{
  int d = deg(p);
  vector<ZZ> v(d+1);
  for(int i=0; i<=d; i++)
    v[i] = coeff(p, i);
  return v;
}

vector<ZZ> coeffs(const ZZ_pX& p)
{
  int d = deg(p);
  vector<ZZ> v(d+1);
  for(int i=0; i<=d; i++)
    v[i] = rep(coeff(p, i));
  return v;
}

string polynomial_string(const ZZX& p, const string& var)
{
  return polynomial_string(coeffs(p), var);
}

string polynomial_string(const ZZ_pX& p, const string& var)
{
  return polynomial_string(coeffs(p), var);
}

template<class T>
mat_ZZ mat_to_mat_ZZ(Zmat<T> A)
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
template mat_ZZ mat_to_mat_ZZ<int>(Zmat<int> A);
template mat_ZZ mat_to_mat_ZZ<long>(Zmat<long> A);
template mat_ZZ mat_to_mat_ZZ<bigint>(Zmat<bigint> A);

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

// Coprime test
int AreCoprime(const ZZX& f, const ZZX& g)
{
  return deg(GCD(f, g)) == 0;
}

// Squarefree test
int IsSquareFree(const ZZX& f)
{
  return AreCoprime(f, diff(f));
}

// Irreducibility test (ignoring content)
int IsIrreducible(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  return factors.length()==1;
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
template mat_m evaluate<bigint>(const ZZX& f, const Zmat<bigint>& A);

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

// display factors of a polynomaial:
void display_factor(const pair_ZZX_long& f)
{
  ZZX p = f.a;
  string pol = polynomial_string(p);
  int d = deg(p), e = f.b;
  cout << "(degree " << d << ")\t"
       << pol
       << "\t to power " << e;
  //cout << " (coefficients " << p << ")";
}

void display_factors(const ZZX& f)
{
  ZZ cont; vec_pair_ZZX_long factors;
  factor(cont, factors, f);
  ::sort(factors.begin(), factors.end(), fact_cmp);
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout << (i+1) << ":\t";
      display_factor(factors[i]);
      cout<<endl;
    }
}

void display_factor(const pair_ZZ_pX_long& f)
{
  ZZ_pX p = f.a;
  string pol = polynomial_string(p);
  int d = deg(p), e = f.b;
  cout << "(degree " << d << ")\t"
       << pol
       << "\t to power " << e;
  //cout << " (coefficients " << p << ")";
}

void display_factors(const ZZ_pX& f)
{
  vec_pair_ZZ_pX_long factors = berlekamp(f);
  ::sort(factors.begin(), factors.end(), fact_modp_cmp);
  long nf = factors.length();
  for(int i=0; i<nf; i++)
    {
      cout << (i+1) << ":\t";
      display_factor(factors[i]);
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
