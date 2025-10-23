// File HECKEFIELD.CC: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#include <NTL/mat_ZZ.h>
#include <NTL/mat_poly_ZZ.h>
#include <NTL/ZZXFactoring.h>
#include "eclib/subspace.h"
#include "matprocs.h"
#include "heckefield.h"

HeckeField::HeckeField(const ZZX& p, string a, int verb)
{
  if (IsMonic(p) && IsIrreducible(p))
    {
      // Set A to be the companion matrix of p.  The function
      // CompanionMatrix(p) returns a mat_ZZ so we do this manually.
      d = deg(p);
      mat_m m(d,d);
      ZZ one(1);
      for(int i=1; i<d; i++)
        {
          m(i+1,i) = one;
          m(i,d) = -coeff(p, i-1);
        }
      m(d,d) = -coeff(p, d-1);
      // Finally call the other constructor
      *this = HeckeField(m, one, a, verb);
    }
}

HeckeField::HeckeField() // defaults to Q
{
  d=1;
  ZZ one(1);
  denom=one;
  A.init(1,1); // zero matrix
  B.init(1,1); // zero matrix
  C.init(1,1); // zero matrix
  B(1,1) = Bdet = Bdet1 = Bdet2 = Bdet3 = one;
  Binv = B;
  var = ""; // will not be used anyway
}

HeckeField::HeckeField(const mat_m& m, const ZZ& den, string a, int verb)
  : var(a), d(m.nrows()), denom(den), A(m)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In HeckeField constructor" << endl;
    }
  minpoly = scaled_charpoly(mat_to_mat_ZZ(A), denom);
  if (verb)
    cout << " - min poly = " << polynomial_string(minpoly) << ", generator " << var << endl;


  // Compute change of basis matrix B, with column j equal to
  // denom^(n-j)*A^(j-1)v for j from 1 to d
  vec_m v(d);
  v[1] = pow(denom,d-1); // so v=[1,0,...,0]*denom^(d-1)
  // cout<<"v = "<<v<<endl;
  B.init(d,d);
  B.setcol(1,v);
  for(int i=2; i<=d; i++)
    {
      v = A*v / denom;
      B.setcol(i,v);
    }
  Bcontent = B.content();
  // cout << "Content of original B is " << Bcontent << endl;
  B /= Bcontent;
  Binv.init(d,d);
  Bdet = inverse(B,Binv); // so B*Binv = Bdet*identity
  Bdet1 = B(1,1);
  Bdet2 = Binv(1,1);
  Bdet3 = Bdet2*denom;
  mat_m I = mat_m::identity_matrix(d);
  assert (B*Binv == Bdet*I);
  // Now we should have Binv*A*B = denom*Bdet * C, where C = companion
  // matrix of minpoly. i.e. the B-conjugate of A can be divided by
  // denom to give the integral matrix C.
  C = (Binv*A*B) / (denom*Bdet);
  assert (CompanionMatrix(minpoly) == mat_to_mat_ZZ(C));

  Cpowers.resize(d);
  Cpowers[0] = I;
  for (int i=1; i<d; i++)
    Cpowers[i] = C*Cpowers[i-1];
  if(verb)
    {
      cout<<"basis  matrix = ";
      output_flat_matrix(Binv);
      cout<<endl;
      cout<<"inverse basis matrix = ";
      output_flat_matrix(B);
      cout<<endl;
      cout << "basis factor  = " << Bdet1 << "*" << Bdet2 << "=" << Bdet <<endl;
      if (verb>1)
        {
          cout<<"companion matrix  = ";
          output_flat_matrix(C);
          cout<<endl;
        }
      cout << "Leaving HeckeField constructor" << endl;
      cout << "----------------------------"<<endl;
    }
}

void HeckeField::display_bases(ostream&s) const
{
  mat_m I(mat_m::identity_matrix(d));

  s << "Powers of A (i.e. powers of " << var << " in A-embedding):\n";
  vector<mat_m> Apowers(d);
  Apowers[0] = I;
  for (int i=1; i<d; i++)
    Apowers[i] = A*Apowers[i-1];
  for (auto Apow: Apowers)
    {
      Apow.output(s);
      s<<endl;
      s<<endl;
    }
  ZZ fac = pow(denom, d-1) * Bdet;
  s << "Basis in A-embedding, scaled by "<< denom << "^" << d-1 << " * " << Bdet << " = " << fac <<":\n";
  fac = Bdet*Bcontent;
  s << "(first columns should be standard basis vectors * "<< fac <<")\n";
  for(int i=1; i<=d; i++)
    {
      vec_m coli = Binv.col(i);
      for(int j=1; j<=d; j++)
        coli[j] *= pow(denom, d-j);
      mat_m M = lin_comb_mats(coli, Apowers);
      M.output(s);
      s<<endl;
      s<<endl;
      if (!(M.col(1)==fac*I.col(i)))
        {
          cout<<"Scaling problem:  column(1) = "<<M.col(1)
              <<" but fac*e_"<<i<<" = "<<fac*I.col(i)<<endl;
        }
    }
  s << "Powers of C (i.e. powers of " << var << " in C-embedding):\n";
  for (auto Cpow: Cpowers)
    {
      Cpow.output(s);
      s<<endl;
      s<<endl;
    }
  s << "Basis in C-embedding, scaled by "<<Bdet<<":\n";
  s << "(first columns should be columns of basis())\n";
  for(int i=1; i<=d; i++)
    {
      mat_m M = lin_comb_mats(Binv.col(i), Cpowers);
      M.output(s);
      s<<endl;
      s<<endl;
      assert (M.col(1) == Binv.col(i));
    }
}

void HeckeField::display(ostream&s, int raw)
{
  static const ZZ one(1);
  string fpol = polynomial_string(minpoly);
  if (d==1)
    {
      s << "Q" << endl;
      return;
    }
  s << "Q(" << var << ") with defining polynomial "<< fpol <<" of degree "<<d;
  if (d==2)
    s << ", discriminant "<<discriminant(minpoly);
  s << endl;
  if(raw)
    {
      s << "   Raw basis with respect to alpha-power basis:\n";
      for(int i=1; i<=d; i++)
        {
          HeckeFieldElement bi = element(vec_m::unit_vector(d,i), one, 1);
          s << "   #"<<i<<": "<< bi << endl;
          s << "   (with char poly = " << polynomial_string(bi.charpoly())
            << " and minpoly = " << polynomial_string(bi.minpoly()) << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////

// raw means the given coords are w.r.t. the B-basis
HeckeFieldElement::HeckeFieldElement( HeckeField* HF, const vec_m& c, const ZZ& d, int raw)
    :F(HF), coords(c), denom(d)
{
  // cout<<"Constructing a HeckeFieldElement in Q("<<F->var<<")\n";
  if (raw)
    {
      coords = F->Binv *coords;
      denom *= F->Bdet3;
    }
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

// creation from a rational
HeckeFieldElement::HeckeFieldElement( HeckeField* HF, const ZZ& a, const ZZ& d)
    :F(HF), coords(a*vec_m::unit_vector(HF->d, 1)), denom(d)
{
  // cout<<"Constructing a HeckeFieldElement in Q("<<F->var<<")\n";
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

HeckeFieldElement HeckeField::element(const vec_m& c, const ZZ& d, int raw)
{
  return HeckeFieldElement(this, c, d, raw);
}

HeckeFieldElement HeckeField::zero()
{
  return HeckeFieldElement(this, vec_m(d));
}

HeckeFieldElement HeckeField::one()
{
  return HeckeFieldElement(this, vec_m::unit_vector(d, 1));
}

HeckeFieldElement HeckeField::gen()
{
  return HeckeFieldElement(this, vec_m::unit_vector(d, 2));
}

void HeckeFieldElement::cancel() // divides through by gcd(content(coords, denom))
{
  if (IsOne(denom))
    return;
  ZZ g = gcd(content(coords), denom);
  if (IsOne(g))
    return;
  denom /=g;
  coords /= g;
}

int HeckeFieldElement::is_zero() const
{
  return trivial(coords);
}

int HeckeFieldElement::is_one() const
{
  return IsOne(denom) && coords == vec_m::unit_vector(F->d,1);
}

int HeckeFieldElement::is_minus_one() const
{
  return IsOne(denom) && coords == -vec_m::unit_vector(F->d,1);
}

string HeckeFieldElement::str() const
{
  // cout << "In HeckeFieldElement::str() with var = " << F->var << endl;
  string s1 = polynomial_string(coords, F->var);
  if (denom==1)
    return s1;
  ostringstream s;
  s << "(" << s1 << ")/" << denom;
  return s.str();
}

int HeckeFieldElement::operator==(const HeckeFieldElement& b) const
{
  return F==b.F && denom==b.denom && coords==b.coords;
}

mat_m HeckeFieldElement::matrix() const // ignores denom
{
  return lin_comb_mats(coords, F->Cpowers);
}

// the charpoly is a power of the minpoly
ZZX HeckeFieldElement::minpoly() const
{
  return factor(charpoly())[0].a;
}

HeckeFieldElement HeckeFieldElement::operator+(const HeckeFieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  return HeckeFieldElement(F, b.denom*coords + denom*b.coords, denom*b.denom);
}

HeckeFieldElement HeckeFieldElement::operator-() const
{
  return HeckeFieldElement(F, -coords, denom);
}

HeckeFieldElement HeckeFieldElement::operator-(const HeckeFieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  return HeckeFieldElement(F, b.denom*coords - denom*b.coords, denom*b.denom);
}

HeckeFieldElement HeckeFieldElement::operator*(const HeckeFieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      exit(1);
    }
  return HeckeFieldElement(F, (matrix()*b.matrix()).col(1), denom*b.denom);
}

HeckeFieldElement HeckeFieldElement::inverse() const // raise error if zero
{
  if (is_zero())
    {
      cerr << "Attempt to invert zero!" << endl;
      exit(1);
    }
  mat_m M = matrix(), Minv;
  ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
  HeckeFieldElement ans(F, denom*Minv.col(1), Mdet);
  assert (operator*(ans).is_one());
  return ans;
}

HeckeFieldElement HeckeFieldElement::operator/(const HeckeFieldElement& b) const // raise error if b is zero
{
  if (b.is_zero())
    {
      cerr << "Attempt to divide by zero!" << endl;
      exit(1);
    }
  if (F!=b.F)
    {
      cerr << "Attempt to divide elements of different fields!" << endl;
      exit(1);
    }
  return operator*(b.inverse());
}
