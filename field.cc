// File FIELD.CC: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#include "field.h"

Field::Field(const ZZX& p, string a, int verb)
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
      *this = Field(m, one, a, verb);
    }
}

Field::Field() // defaults to Q
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

Field::Field(const mat_m& m, const ZZ& den, string a, int verb)
  : var(a), d(m.nrows()), denom(den), A(m)
{
  if (verb)
    {
      cout << "----------------------------"<<endl;
      cout << "In Field constructor" << endl;
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
      cout << "Leaving Field constructor" << endl;
      cout << "----------------------------"<<endl;
    }
}

void Field::display_bases(ostream&s) const
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

void Field::display(ostream&s, int raw)
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
          FieldElement bi = element(vec_m::unit_vector(d,i), one, 1);
          s << "   #"<<i<<": "<< bi << endl;
          s << "   (with char poly = " << polynomial_string(bi.charpoly())
            << " and minpoly = " << polynomial_string(bi.minpoly()) << "\n";
        }
    }
}

////////////////////////////////////////////////////////////////////////

// raw means the given coords are w.r.t. the B-basis
FieldElement::FieldElement( Field* HF, const vec_m& c, const ZZ& d, int raw)
    :F(HF), coords(c), denom(d)
{
  // cout<<"Constructing a FieldElement in Q("<<F->var<<")\n";
  if (raw)
    {
      coords = F->Binv *coords;
      denom *= F->Bdet3;
    }
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

// creation from a rational
FieldElement::FieldElement( Field* HF, const ZZ& a, const ZZ& d)
    :F(HF), coords(a*vec_m::unit_vector(HF->d, 1)), denom(d)
{
  // cout<<"Constructing a FieldElement in Q("<<F->var<<")\n";
  cancel();
  // cout<<" - finished constructing " << (*this) << "\n";
}

FieldElement Field::element(const vec_m& c, const ZZ& d, int raw)
{
  return FieldElement(this, c, d, raw);
}

FieldElement Field::zero()
{
  return FieldElement(this, vec_m(d));
}

FieldElement Field::one()
{
  return FieldElement(this, vec_m::unit_vector(d, 1));
}

FieldElement Field::gen()
{
  return FieldElement(this, vec_m::unit_vector(d, 2));
}

void FieldElement::cancel() // divides through by gcd(content(coords, denom))
{
  if (IsOne(denom))
    return;
  ZZ g = gcd(content(coords), denom);
  if (IsOne(g))
    return;
  denom /=g;
  coords /= g;
}

int FieldElement::is_zero() const
{
  return trivial(coords);
}

int FieldElement::is_one() const
{
  return IsOne(denom) && coords == vec_m::unit_vector(F->d,1);
}

int FieldElement::is_minus_one() const
{
  return IsOne(denom) && coords == -vec_m::unit_vector(F->d,1);
}

string FieldElement::str() const
{
  // cout << "In FieldElement::str() with var = " << F->var << endl;
  string s1 = polynomial_string(coords, F->var);
  if (denom==1)
    return s1;
  ostringstream s;
  s << "(" << s1 << ")/" << denom;
  return s.str();
}

int FieldElement::operator==(const FieldElement& b) const
{
  return F==b.F && denom==b.denom && coords==b.coords;
}

mat_m FieldElement::matrix() const // ignores denom
{
  return lin_comb_mats(coords, F->Cpowers);
}

// the charpoly is a power of the minpoly
ZZX FieldElement::minpoly() const
{
  return factor(charpoly())[0].a;
}

// add b to this
void FieldElement::operator+=(const FieldElement& b)
{
  if (F!=b.F)
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  coords = b.denom*coords + denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator+(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to add elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  a += b;
  return a;
}

FieldElement FieldElement::operator-() const
{
  return FieldElement(F, -coords, denom);
}

// subtract b
void FieldElement::operator-=(const FieldElement& b)
{
  if (F!=b.F)
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  coords = b.denom*coords - denom*b.coords;
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator-(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to subtract elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  a -= b;
  return a;
}

void FieldElement::operator*=(const FieldElement& b) // multiply by b
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      exit(1);
    }
  coords = (matrix()*b.matrix()).col(1);
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator*(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      exit(1);
    }
  FieldElement a = *this;
  a *= b;
  return a;
}

FieldElement FieldElement::inverse() const // raise error if zero
{
  if (is_zero())
    {
      cerr << "Attempt to invert zero!" << endl;
      exit(1);
    }
  mat_m M = matrix(), Minv;
  ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
  FieldElement ans(F, denom*Minv.col(1), Mdet);
  assert (operator*(ans).is_one());
  return ans;
}

void FieldElement::operator/=(const FieldElement& b)      // divide by b
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
  operator*=(b.inverse());
}

FieldElement FieldElement::operator/(const FieldElement& b) const // raise error if b is zero
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
  FieldElement a = *this;
  a /= b;
  return a;
}

FieldElement evaluate(const ZZX& f, const FieldElement a)
{
  long d = deg(f);
  FieldElement fa = a;
  fa *= coeff(f,d);
  for(int i=d-1; i>=0; i--)
    {
      fa += coeff(f,i);
      if(i)
        fa *= a;
    }
  return fa;

}

// NB for a in F, either [Q(sqrt(a))=Q(a)] or [Q(sqrt(a)):Q(a)]=2.
// The first function only applies when a has maximal degree:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
int FieldElement::is_absolute_square(FieldElement& r) const
{
  return is_absolute_square(r, minpoly());
}

// Same as above if the min poly is known
int FieldElement::is_absolute_square(FieldElement& r, const ZZX& minpol)  const
{
  ZZX g, g0, g1;
  if (::is_square(minpol, g))
    {
      parity_split(g, g0, g1);
      // Now 0 = g(r) = g0(a)+r*g1(a) and g1(a)!=0
      r = - evaluate(g0,*this)/evaluate(g1,*this);
      if (r*r== *this)
        return 1;
      else
        {
          cout << (*this) << " has min poly " << polynomial_string(minpol) << " whose double has factor "
               << polynomial_string(g)
               << " with even part g0 = " << polynomial_string(g0)
               << " and odd part g1 = " << polynomial_string(g1) << endl;
          cout << "These evaluate to " << evaluate(g0,*this) << " and " << evaluate(g1,*this)
               << " with negative quotient r = " << r
               << " but r*r = " << r*r << endl;
          return 0;
        }
    }
  return 0;
}

// The second function applies in general:
// return 1 and r s.t. r^2=this, with deg(r)=degree(), else 0
int FieldElement::is_square(FieldElement& r, int ntries) const
{
  if (is_zero() || is_one())
    {
      r = *this;
      return 1;
    }
  // If this has full degree, or if the relative degree is odd, just
  // call is_absolute_square()
  ZZX m = minpoly();
  if (((F->d)/deg(m))%2)
    {
      //cout << "Is_square() succeeds directly" << endl;
      return is_absolute_square(r, m);
    }

  // Otherwise we multiply this by a "random" square to have full degree first
  FieldElement b = F->gen();
  for (int i=0; i<ntries; i++, b+=to_ZZ(1))
    {
      FieldElement abb = (*this)*b*b, rb(F);
      ZZX mb = abb.minpoly();
      if (((F->d)/deg(mb))%2)
        {
          if (is_absolute_square(rb, mb))
            {
              //cout << "Is_square() succeeds after " << i+1 << " tries" << endl;
              r = rb/b;
              return 1;
            }
          else
            return 0;
        }
      // else keep trying
    } // end of loop over shifts i
  cout << "is_square() fails on " << (*this) << " after " << ntries << " tries" << endl;
  return 0;
}
