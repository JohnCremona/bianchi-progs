// File FIELD.CC: class for working with number fields for Hecke eigenvalues
//////////////////////////////////////////////////////////////////////////

#include "field.h"

Field* FieldQQ = new Field();

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
  if (denom<0)
    {
      denom = -denom;
      coords = -coords;
    }
  ZZ g = gcd(content(coords), denom);
  if (IsOne(g))
    return;
  denom /=g;
  coords /= g;
}

int FieldElement::is_zero() const
{
  if (F->isQ())
    return ::is_zero(coords[1]);
  return trivial(coords);
}

int FieldElement::is_one() const
{
  if (F->isQ())
    return IsOne(coords[1]) && IsOne(denom);
  return IsOne(denom) && coords == vec_m::unit_vector(F->d,1);
}

int FieldElement::is_minus_one() const
{
  if (F->isQ())
    return IsOne(-coords[1]) && IsOne(denom);
  return IsOne(denom) && coords == -vec_m::unit_vector(F->d,1);
}

string FieldElement::str() const
{
  // cout << "In FieldElement::str() with var = " << F->var << endl;
  // cout << "coords = " << coords << ", denom = " << denom << endl;
  ostringstream s;
  if (F->isQ())
    {
      s << coords[1];
      if (denom!=1)
        s << "/" << denom;
      return s.str();
    }
  string n = polynomial_string(coords, F->var);
  if (denom==1)
    return n;
  s << "(" << n << ")/" << denom;
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

// the charpoly is a power of the irreducible minpoly NB This monic
// integer polynomial is the min poly of the numerator, it ignoes the
// denominator!
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
  if (b.is_zero())
    return;
  if (F->isQ())
    coords[1] = b.denom*coords[1] + denom*b.coords[1];
  else
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
  if (!b.is_zero())
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
  if (b.is_zero())
    return;
  if (F->isQ())
    coords[1] = b.denom*coords[1] - denom*b.coords[1];
  else
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
  if (!b.is_zero())
    a -= b;
  return a;
}

void FieldElement::operator*=(const FieldElement& b) // multiply by b
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (is_zero()) return;
  if (b.is_zero()) {*this = b; return;}
  if (F->isQ())
    coords[1] *= b.coords[1];
  else
    coords = (matrix()*b.matrix()).col(1);
  denom *= b.denom;
  cancel();
}

FieldElement FieldElement::operator*(const FieldElement& b) const
{
  if (F!=b.F)
    {
      cerr << "Attempt to multiply elements of different fields!" << endl;
      cerr << "LHS field: "; F->display(); cerr << endl;
      cerr << "RHS field: "; b.F->display(); cerr << endl;
      exit(1);
    }
  if (b.is_zero()) return b;
  FieldElement a = *this;
  if (!is_zero())
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
  if (is_one() || is_minus_one())
    {
      return FieldElement(*this);
    }
  FieldElement ans(F);
  if (F->isQ())
    {
      ans = FieldElement(F, denom, coords[1]);
      // constructor will ensure new denom>0
    }
  else
    {
      mat_m M = matrix(), Minv;
      ZZ Mdet = ::inverse(M,Minv); // so M*Minv = Mdet*identity
      ans = FieldElement(F, denom*Minv.col(1), Mdet);
    }
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
  if (is_zero())
    return;
  if (F->isQ())
    {
      coords[1] *= b.denom;
      denom *= b.coords[1];
      cancel();
    }
  else
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
  if (!is_zero())
    a /= b;
  return a;
}

FieldElement evaluate(const ZZX& f, const FieldElement a)
{
  FieldElement fa(a.F, to_ZZ(0));
  for(int i=deg(f); i>=0; i--)
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
  if (F->isQ())
    {
      ZZ n = coords[1];
      if (n<0)
        return 0;
      ZZ x = n*denom, rn;
      SqrRoot(rn, x); // rounds down if x is not square;
      if (rn*rn != x)
          return 0;
      // sqrt(n/d) = sqrt(n*d)/d = r/d
      r = FieldElement(F, rn, denom);
      assert (r*r== *this);
      return 1;
    }
  return is_absolute_square(r, minpoly());
}

// Same as above if the min poly is known
int FieldElement::is_absolute_square(FieldElement& r, const ZZX& minpol)  const
{
  ZZX g, g0, g1;
  ZZX f = scale_poly_up(minpol, denom);
  if (::is_square(f, g))
    {
      parity_split(g, g0, g1);
      // Now 0 = g(r) = g0(a)+r*g1(a) and g1(a)!=0
      r = - evaluate(g0,*this)/(evaluate(g1,*this)*denom);
      if (r*r== *this)
        return 1;
      else
        {
          cout << (*this) << " has scaled min poly " << polynomial_string(f)
               << " whose double has factor " << polynomial_string(g)
               << " with even part g0 = " << polynomial_string(g0)
               << " and odd part g1 = " << polynomial_string(g1) << endl;
          cout << "These evaluate to " << evaluate(g0,*this) << " and " << evaluate(g1,*this)
               << " with negative quotient/denom r = " << r
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
  if (F->isQ())
    {
      return is_absolute_square(r);
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

unsigned int FieldModSq::get_index(const FieldElement& a, FieldElement& s)
{
  unsigned int i=0;
  for (auto x: elements)
    {
      if ((a*x).is_square(s))
        {
          assert (a*elements[i] == s*s);
          // Now a*x = s^2 so a = x*(s/x)^2
          // but we want a = x*s^2
          s /= x;
          assert (a == elements[i]*s*s);
          return i;
        }
      i++;
    }
  // We get here if a is not in the current group. For fields other
  // than Q we use a itself as the new gen, otherwise the squarefree
  // part of a.
  FieldElement newgen = a;
  ZZ a1;
  if (F->isQ())
    {
      a1 = squarefree_part(a.coords[1]*a.denom);
      newgen = FieldElement(F, a1);
    }
  gens.push_back(newgen);
  r++;
  s = F->one();
  // cout << "rank of k*/(k*)^2 grows to " << r << " after adding generator " << newgen << endl;
  // cout << "elements were " << elements << endl;
  vector<FieldElement> new_elements(elements.size(), FieldElement(F));
  if (F->isQ())
    {
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [this, a1](const FieldElement& x){return FieldElement(F, squarefree_part(a1*x.coords[1]*x.denom));});
    }
  else
    {
      std::transform(elements.begin(), elements.end(), new_elements.begin(),
                     [newgen](const FieldElement& x){return newgen*x;});
    }
  elements.insert(elements.end(), new_elements.begin(), new_elements.end());
  // cout << "elements are now " << elements << endl;
  return r;
}

string FieldModSq::elt_str(unsigned int i) const
{
  static const string v = "r";
  ostringstream s;
  if (i==0)
    {
      s << "1";
    }
  else
    {
      int first = 1;
      for (unsigned int j=0; j<r; j++)
        {
          if (bit(i,j))
            {
              if (!first) s << "*";
              s << v << (j+1);
              first = 0;
            }
        }
    }
  return s.str();
}

void FieldModSq::display()
{
  cout << "FieldModSq data:\n";
  cout << "Field  is " <<flush;
  F->display();
  cout << "rank = " << r;
  cout << ", gens = " << gens << endl;
  cout << "order = " << elements.size();
  cout << ", elements = " << elements << endl;
}

Eigenvalue Eigenvalue::operator*(Eigenvalue b) const
{
  if (is_zero()) return Eigenvalue(*this);
  if (b.is_zero()) return b;
  FieldElement r = root_part() * b.root_part();
  FieldElement s(a.F);
  unsigned int j = SqCl->get_index(r, s);
  assert (r == s*s*SqCl->elt(j));
  // Now r = s^2 * elt(j)
  return Eigenvalue(a*b.a*s, SqCl, j);
}

Eigenvalue Eigenvalue::operator/(Eigenvalue b) const
{
  if (is_zero()) return Eigenvalue(*this);
  // cout << "\nDividing " << (*this) << " by " << b << endl;
  FieldElement r = root_part() / b.root_part();
  // cout << "Quotent of root_parts " << root_part() << " and " << b.root_part() << " is " << r << endl;
  FieldElement s(a.F); // value ignored, just to set the field
  unsigned int j = SqCl->get_index(r, s);
  assert (r == s*s*SqCl->elt(j));
  // cout << "... which has index " << j << ", elt(j) = "<< SqCl->elt(j)
  //      << " with s="<<s<<endl;
  // cout << "s^2*elt(j) = " << s*s*SqCl->elt(j) << endl;
  // Now r = s^2 * elt(j)
  return Eigenvalue(s*a/b.a, SqCl, j);
}

string Eigenvalue::str() const
{
  ostringstream s;
  if (is_zero())
    s << "0";
  else
    {
      int QQ = a.F->isQ();
      if (!root_index) // then a involves no sqrts
        s << a;
      else // it does
        {
          if (QQ)
            {
              if (a.is_one())
                ;
              else
                if (a.is_minus_one())
                  s << "-";
                else
                  s << a << "*";
            }
          else
            {
              s << "(" << a << ")";
            }
          s << "sqrt(";
          if (QQ)
            s << SqCl->elt(root_index);
          else
            s << SqCl->elt_str(root_index);
          s  << ")";
        }
    }
  return s.str();
}
